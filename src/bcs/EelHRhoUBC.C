#include "EelHRhoUBC.h"

template<>
InputParameters validParams<EelHRhoUBC>()
{
  InputParameters params = validParams<IntegratedBC>();
    params.addRequiredParam<std::string>("equation_name", "The name of the equation this BC is acting on");
    params.addParam<bool>("isLiquid", true, "boolean for phase");
    // Coupled conservative variables:
    params.addCoupledVar("alA", "liquid volume fraction: alphaA");
    params.addCoupledVar("alrhoA", "density: rhoA");
    params.addCoupledVar("alrhouA_x", "x component of the momentum: rhouA_x");
    params.addCoupledVar("alrhoEA", "total energy: rho*E*A");
    // Coupled aux variables:
    params.addCoupledVar("area", "Coupled area variable");
    // Input parameters
    params.addRequiredParam<Real>("rhou_liquid", "Specified momentum for liquid: rhou");
    params.addRequiredParam<Real>("H_liquid", "Specified enthalpy for liquid");
    params.addRequiredParam<Real>("rhou_vapor", "Specified momentum for vapor: rhou");
    params.addRequiredParam<Real>("H_vapor", "Specified enthalpy for vapor");
    params.addRequiredParam<Real>("alpha_bc", "Specified liquid volume fraction");
    // Make the name of the EOS function a required parameter.
    params.addRequiredParam<UserObjectName>("eos", "The name of equation of state object to use.");

  return params;
}

EelHRhoUBC::EelHRhoUBC(const std::string & name, InputParameters parameters) :
    IntegratedBC(name, parameters),
    // Type of equation:
    _eqn_name(getParam<std::string>("equation_name")),
    _eqn_type("VOIDFRACTION, CONTINUITY, XMOMENTUM, ENERGY, INVALID", _eqn_name),
    _isLiquid(getParam<bool>("isLiquid")),
    // Coupled aux variables:
    _alA(coupledValue("alA")),
    _alrhoA(coupledValue("alrhoA")),
    _alrhouA_x(coupledValue("alrhouA_x")),
    _alrhoEA(coupledValue("alrhoEA")),
    // Coupled aux variables:
    _area(coupledValue("area")),
    // Stagnation variables:
    _rhou(_isLiquid ? getParam<Real>("rhou_liquid") : getParam<Real>("rhou_vapor")),
    _H(_isLiquid ? getParam<Real>("H_liquid") : getParam<Real>("H_vapor")),
    _alpha_bc(_isLiquid ? getParam<Real>("alpha_bc") : 1.-getParam<Real>("alpha_bc")),
    // Equation of state:
    _eos(getUserObject<EquationOfState>("eos"))
    // Parameters for jacobian matrix:
//    _rhoA_nb(coupled("rhoA")),
//    _rhouA_x_nb(coupled("rhouA_x")),
//    _rhouA_y_nb(isCoupled("rhouA_y") ? coupled("rhouA_x") : -1),
//    _rhoEA_nb(coupled("rhoEA"))
{}

Real
EelHRhoUBC::computeQpResidual()
{
  Real rho, P;
  switch (_eqn_type)
  {
  case CONTINUITY:
        return _alpha_bc * _area[_qp] * _rhou * _normals[_qp](0) * _test[_i][_qp];
  case XMOMENTUM:
        rho = _alrhoA[_qp] / (_alpha_bc*_area[_qp]);
          P = _eos.pressure(rho, _alrhouA_x[_qp]/_alrhoA[_qp], _alrhoEA[_qp]/(_alpha_bc*_area[_qp]));
        return _alpha_bc * _area[_qp] * ( _rhou * _rhou / rho + P ) * _normals[_qp](0) * _test[_i][_qp];
  case ENERGY:
        return _alpha_bc * _area[_qp] * _rhou * _H * _normals[_qp](0) * _test[_i][_qp];
  default:
    mooseError("The equation with name: \"" << _eqn_name << "\" is not supported in the \"EelHRhoUBC\" type of boundary condition.");
  }
}

Real
EelHRhoUBC::computeQpJacobian()
{
    return 0.;
}

Real
EelHRhoUBC::computeQpOffDiagJacobian(unsigned _jvar)
{
    return 0.;
}
