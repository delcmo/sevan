#include "EelStagnationPandTBC.h"

template<>
InputParameters validParams<EelStagnationPandTBC>()
{
  InputParameters params = validParams<IntegratedBC>();

  params.addRequiredParam<std::string>("equation_name", "The name of the equation this BC is acting on");

    // Coupled conservative variables:
    params.addCoupledVar("alrhoA", "alpha*rho*A");
    params.addCoupledVar("alrhouA_n", "component of the momentum normal to the surface: alpha*rho*u*A");
    // Coupled aux variables:
    params.addCoupledVar("area", "Coupled area variable");
    // Input parameters
    params.addRequiredParam<Real>("p0_bc", "Stagnation pressure at the boundary");
    params.addRequiredParam<Real>("T0_bc", "Liquid stagnation temperature at the boundary");
    params.addParam<Real>("gamma0_bc", 0., "Stagnation angle");
    params.addParam<Real>("alpha0_bc", 1., "Liquid void fraction at the boundary");
    // Make the name of the EOS function a required parameter.
    params.addRequiredParam<UserObjectName>("eos", "The name of equation of state object to use.");
    // Boolean
    params.addParam<bool>("isLiquid", true, "is liquid or not?");
    params.addParam<bool>("is5EquModel", false, "is 5 equations model?");

  return params;
}

EelStagnationPandTBC::EelStagnationPandTBC(const std::string & name, InputParameters parameters) :
    IntegratedBC(name, parameters),
    // Type of equation:
    _eqn_name(getParam<std::string>("equation_name")),
    _eqn_type("VOIDFRACTION, CONTINUITY, XMOMENTUM, YMOMENTUM, ZMOMENTUM, ENERGY, VOID_FRACTION, INVALID", _eqn_name),
    // Coupled aux variables:
    _alrhoA(coupledValue("alrhoA")),
    _alrhouA_n(coupledValue("alrhouA_n")),
    // Coupled aux variables:
    _area(coupledValue("area")),
    // Stagnation variables:
    _p0_bc(getParam<Real>("p0_bc")),
    _T0_bc(getParam<Real>("T0_bc")),
    _gamma0_bc(getParam<Real>("gamma0_bc")),
    _alpha_bc_l(getParam<Real>("alpha0_bc")),
    // Equation of state:
    _eos(getUserObject<EquationOfState>("eos")),
    // Boolean:
    _isLiquid(getParam<bool>("isLiquid"))
{
    // Pre-compute some stagnation cpefficients:
    _rho0_bc = _eos.rho_from_p_T(_p0_bc, _T0_bc);
    _H0_bc = _eos.e_from_p_rho(_p0_bc, _rho0_bc) + _p0_bc / _rho0_bc;
    _K = (_p0_bc + _eos.Pinf()) / std::pow(_rho0_bc, _eos.gamma());
    _H_bar = _eos.gamma() * (_p0_bc + _eos.Pinf()) / _rho0_bc / (_eos.gamma() - 1);
}

Real
EelStagnationPandTBC::computeQpResidual()
{
    // Compute the void fraction:
    Real alpha_bc = _isLiquid ? _alpha_bc_l : (1-_alpha_bc_l);
    
    // Compute u_star and v_star:
    Real u_star = _alrhouA_n[_qp]/_alrhoA[_qp];
    Real v_star = u_star * std::tan(_gamma0_bc);
    Real _norm_vel_star2 = u_star*u_star + v_star*v_star;
    
    // Compute rho_star and static pressure:
    Real rho_star = std::pow((_H_bar - 0.5*_norm_vel_star2)*(_eos.gamma()-1)/(_eos.gamma())/_K, 1./(_eos.gamma()-1));
    Real p_bc = _K * std::pow(rho_star, _eos.gamma()) - _eos.Pinf();
    
  switch (_eqn_type)
  {
      case CONTINUITY:
//          return _alrhouA_n[_qp] * _normals[_qp](0) * _test[_i][_qp];
          return alpha_bc * rho_star * u_star * _area[_qp] * _normals[_qp](0) * _test[_i][_qp];
//        return _alpha*_area[_qp]*rho_star*( u_star*_normals[_qp](0) + v_star*_normals[_qp](1) ) * _test[_i][_qp];
      case XMOMENTUM:
//          return (_u[_qp] * _u[_qp] / _alrhoA[_qp] + alpha_bc * _area[_qp] * p_bc) * _normals[_qp](0) * _test[_i][_qp];
          return alpha_bc * _area[_qp] * (rho_star * u_star * u_star + p_bc) * _normals[_qp](0) * _test[_i][_qp];
//        return alpha_bc*_area[_qp]*(u_star*rho_star*(u_star+v_star) + p_bc) * _normals[_qp](0) * _test[_i][_qp];
      case YMOMENTUM:
        return alpha_bc*_area[_qp]*(v_star*rho_star*(u_star+v_star) + p_bc) * _normals[_qp](1) * _test[_i][_qp];
      case ENERGY:
//          return _alrhouA_n[_qp] * _H0_bc * _normals[_qp](0) * _test[_i][_qp];
          return alpha_bc * rho_star * u_star * _area[_qp] * _H0_bc * _normals[_qp](0) * _test[_i][_qp];
//        return alpha_bc*_area[_qp]*rho_star*_H0_bc*(u_star*_normals[_qp](0)+v_star*_normals[_qp](1)) * _test[_i][_qp];
      case VOID_FRACTION:
        return 0.;
      default:
        mooseError("The equation with name: \"" << _eqn_name << "\" is not supported in the \"OneD7EqnStagnationPandTBC\" type of boundary condition.");
        return 0.;
  }
}

Real
EelStagnationPandTBC::computeQpJacobian()
{
  // TODO
  return 0;
}

Real
EelStagnationPandTBC::computeQpOffDiagJacobian(unsigned jvar)
{
  // TODO
  return 0;
}
