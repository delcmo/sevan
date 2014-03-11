#include "MassInflow.h"

template<>
InputParameters validParams<MassInflow>()
{
  InputParameters params = validParams<IntegratedBC>();
    params.addRequiredParam<std::string>("equation_name", "The name of the equation this BC is acting on");
    // Coupled variables:
    params.addRequiredCoupledVar("alA", "alpha*A");
    params.addRequiredCoupledVar("alrhoA", "alpha*rho*A");
    params.addRequiredCoupledVar("alrhouA", "alpha*rho*u*A");
    params.addCoupledVar("alrhovA", "alpha*rho*v*A");
    params.addRequiredCoupledVar("alrhoEA", "alpha*rho*E*A");
    // Coupled aux variables:
    params.addCoupledVar("area", "area aux variable");
    // Input parameters:
    params.addRequiredParam<Real>("mass_inflow", "Specified inflow momentum: rho*u*A.");
    params.addRequiredParam<Real>("temp_inflow", "Specified inflow temperature.");
    params.addParam<Real>("angle_inflow", 0., "Inflow angle: only required for 2D runs");
    // Equation of state:
    params.addRequiredParam<UserObjectName>("eos", "The name of equation of state object to use.");
    // Phase
    params.addParam<bool>("isLiquid", true, "boolean");
  return params;
}

MassInflow::MassInflow(const std::string & name, InputParameters parameters) :
    IntegratedBC(name, parameters),
    // Name of the equation:
    _eqn_name(getParam<std::string>("equation_name")),
    _eqn_type("CONTINUITY, XMOMENTUM, YMOMENTUM, ENERGY, INVALID", _eqn_name),
    // Phase:
    _isLiquid(getParam<bool>("isLiquid")),
    // Coupled variables:
    _alA(coupledValue("alA")),
    _alrhoA(coupledValue("alrhoA")),
    _alrhouA(coupledValue("alrhouA")),
    _alrhovA(isCoupled("alrhovA") ? coupledValue("alrhovA") : _zero),
    _alrhoEA(coupledValue("alA")),
    // Coupled aux variables:
    _area(coupledValue("area")),
    // Boundary condition parameters:
    _momA_bc(getParam<Real>("mass_inflow")),
    _T_bc(getParam<Real>("temp_inflow")),
    _theta_bc(getParam<Real>("angle_inflow")),
    // Equation of state:
    _eos(getUserObject<EquationOfState>("eos"))
{
}

Real
MassInflow::computeQpResidual()
{
    // Compute the pressure and void-fraction:
    RealVectorValue vel(_alrhouA[_qp]/_alrhoA[_qp], _alrhovA[_qp]/_alrhoA[_qp], 0.);
    Real press = _eos.pressure(_alrhoA[_qp]/_alA[_qp], vel.size(), _alrhoEA[_qp]/_alA[_qp]);
    Real alpha = _isLiquid ? _alA[_qp] / _area[_qp] : 1. - _alA[_qp] / _area[_qp];
    
    // Compute the boundary values:
    Real rho_bc = _eos.rho_from_p_T(press, _T_bc);
    Real vel_x_bc = cos(_theta_bc)*_momA_bc/(rho_bc*_area[_qp]);
    Real vel_y_bc = sin(_theta_bc)*_momA_bc/(rho_bc*_area[_qp]);
    Real e_bc = _eos.e_from_p_rho(press, rho_bc);
    RealVectorValue vel_bc(vel_x_bc, vel_y_bc, 0.);
    Real E_bc = e_bc + 0.5*vel_bc.size_sq();
    
    std::cout<<"press="<<press<<std::endl;
    std::cout<<"rho="<<rho_bc<<std::endl;
    std::cout<<"vel_x="<<vel_x_bc<<std::endl;
    std::cout<<"vel_y="<<vel_y_bc<<std::endl;
    std::cout<<"mass inflow="<<_momA_bc<<std::endl;
    
    // Switch statement on the equation type:
    switch (_eqn_type) {
        case CONTINUITY:
            return alpha*_momA_bc*_normals[_qp](0)*_test[_i][_qp];
        case XMOMENTUM:
            return _alA[_qp]*( rho_bc*vel_x_bc*vel_bc*_normals[_qp] + press*_normals[_qp](0) )*_test[_i][_qp];
        case YMOMENTUM:
            return _alA[_qp]*( rho_bc*vel_y_bc*vel_bc*_normals[_qp] + press*_normals[_qp](1) )*_test[_i][_qp];
        case ENERGY:
            return _alA[_qp]*( rho_bc*E_bc+press )*vel_bc*_normals[_qp]*_test[_i][_qp];
        default:
            mooseError("The equation with name: \"" << _eqn_name << "\" is not supported in the \"MassInflow\" type of boundary condition.");
        }
}

Real
MassInflow::computeQpJacobian()
{
  // TODO
  return 0;
}

Real
MassInflow::computeQpOffDiagJacobian(unsigned jvar)
{
  // TODO
  return 0;
}
