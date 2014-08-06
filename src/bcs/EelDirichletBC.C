#include "EelDirichletBC.h"

template<>
InputParameters validParams<EelDirichletBC>()
{
  InputParameters params = validParams<NodalBC>();
    params.addRequiredParam<std::string>("equation_name", "The name of the equation this BC is acting on");
    params.addParam<bool>("isLeftBC", true, "boundary: left or right");
    params.addParam<bool>("isLiquid", true, "boolean for phase");
    params.addCoupledVar("area", 1., "Coupled area variable");
    // Initial conditions:
    params.addRequiredParam<Real>("pressure_init_left", "Initial pressure on the left");
    params.addRequiredParam<Real>("pressure_init_right", "Initial pressure on the right");
    params.addRequiredParam<Real>("vel_init_left", "Initial velocity on the left");
    params.addRequiredParam<Real>("vel_init_right", "Inital velocity on the right");
    params.addParam<Real>("temp_init_left", "Initial value of the LIQUID temperature");
    params.addParam<Real>("temp_init_right", "Initial value of the LIQUID temperature");
    params.addParam<Real>("temp_init_left_gas", "Initial value of the GAS temperature if different from liquid temperature");
    params.addParam<Real>("temp_init_right_gas", "Initial value of the GAS temperature if different from liquid temperature");
    params.addParam<Real>("rho_init_left_liq", "Initial value of the LIQUID density");
    params.addParam<Real>("rho_init_right_liq", "Initial value of the LIQUID density");
    params.addParam<Real>("rho_init_left_gas", "Initial value of the GAS density");
    params.addParam<Real>("rho_init_right_gas", "Initial value of the GAS density");
    params.addRequiredParam<Real>("alpha_init_left", "Initial value of the LIQUID void fraction");
    params.addRequiredParam<Real>("alpha_init_right", "Initial value of the LIQUID void fraction");
    // Make the name of the EOS function a required parameter.
    params.addRequiredParam<UserObjectName>("eos", "The name of equation of state object to use.");
  return params;
}

EelDirichletBC::EelDirichletBC(const std::string & name, InputParameters parameters)
  :NodalBC(name, parameters),
    // Type of equation:
    _eqn_type("VOIDFRACTION, CONTINUITY, XMOMENTUM, ENERGY, INVALID", getParam<std::string>("equation_name")),
    // Boundary
    _isLeftBC(getParam<bool>("isLeftBC")),
    _isLiquid(getParam<bool>("isLiquid")),
    // Coupled aux variables:
    _area(coupledValue("area")),
    // Input variables used for ICs:
    _rho(_isLiquid ? (_isLeftBC ? getParam<Real>("rho_init_left_liq") : getParam<Real>("rho_init_right_liq")) : (_isLeftBC ? getParam<Real>("rho_init_left_gas") : getParam<Real>("rho_init_right_gas")) ),
    _temp(_isLiquid ? (_isLeftBC ? getParam<Real>("temp_init_left") : getParam<Real>("temp_init_right")) : (_isLeftBC ? getParam<Real>("temp_init_left_gas") : getParam<Real>("temp_init_right_gas")) ),
    _pressure(_isLeftBC ? getParam<Real>("pressure_init_left") : getParam<Real>("pressure_init_right")),
    _vel(_isLeftBC ? getParam<Real>("vel_init_left") : getParam<Real>("vel_init_right")),
    _alpha(_isLiquid ? (_isLeftBC ? getParam<Real>("alpha_init_left") : getParam<Real>("alpha_init_right")) : (_isLeftBC ? 1.-getParam<Real>("alpha_init_left") : 1.-getParam<Real>("alpha_init_right")) ),
    // Equation of state:
    _eos(getUserObject<EquationOfState>("eos"))
{
    if ( isParamValid("temp_init_left") && isParamValid("temp_init_right") )
        _isDensity = false;
    else
        _isDensity = true;
}

Real
EelDirichletBC::computeQpResidual()
{
    Real density = _rho;
    if (!_isDensity)
        density = (_pressure + _eos.Pinf()) / (_eos.Cv()*(_eos.gamma()-1)*_temp);
    // Total energy
    Real int_energy = (_pressure+_eos.gamma()*_eos.Pinf())/(density*(_eos.gamma()-1)) + _eos.qcoeff();
    Real tot_energy = density*(int_energy + 0.5*_vel*_vel);
    
    switch (_eqn_type)
    {
        case VOIDFRACTION:
            return _u[_qp] - _alpha*_area[_qp];
        case CONTINUITY:
            return _u[_qp] - _alpha*density*_area[_qp];
        case XMOMENTUM:
            return _u[_qp] - _alpha*_area[_qp]*density*_vel;
        case ENERGY:
            return _u[_qp] - _alpha*_area[_qp]*tot_energy;
        default:
            mooseError("The equation with name: \"" << _eqn_type << "\" is not supported in the \"EelDirichletBC\" type of boundary condition.");
    }
}

Real
EelDirichletBC::computeQpJacobian()
  {
    return 0;
    // FIXME: !!!
//    return _phi[_j][_qp];
  }
