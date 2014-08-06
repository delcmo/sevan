#include "EelNeumannBC.h"

template<>
InputParameters validParams<EelNeumannBC>()
{
  InputParameters params = validParams<IntegratedBC>();
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

EelNeumannBC::EelNeumannBC(const std::string & name, InputParameters parameters) :
    IntegratedBC(name, parameters),
    // Type of equation:
    _eqn_name(getParam<std::string>("equation_name")),
    _eqn_type("VOIDFRACTION, CONTINUITY, XMOMENTUM, ENERGY, INVALID", _eqn_name),
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
EelNeumannBC::computeQpResidual()
{
//    std::cout<<"&&&&&&&&&&&&"<<std::endl;
//    _isLeftBC ? std::cout<<"bc=LEFT"<<std::endl : std::cout<<"bc=RIGHT"<<std::endl;
//    _isLiquid ? std::cout<<"LIQDUID" : std::cout<<"VAPOR"<<std::endl;
//    std::cout<<"P="<<_pressure<<std::endl;
//    std::cout<<"v="<<_vel<<std::endl;
//    std::cout<<"rho="<<_rho<<std::endl;
//    std::cout<<"alpha="<<_alpha<<std::endl;
    // Density
    Real density = _rho;
    if (!_isDensity)
        density = (_pressure + _eos.Pinf()) / (_eos.Cv()*(_eos.gamma()-1)*_temp);
    // Total energy
    Real int_energy = (_pressure+_eos.gamma()*_eos.Pinf())/(density*(_eos.gamma()-1)) + _eos.qcoeff();
    Real tot_energy = density*(int_energy + 0.5*_vel*_vel);
    
    switch (_eqn_type)
    {
    case CONTINUITY:
        return _alpha*density*_vel*_area[_qp]*_normals[_qp](0) * _test[_i][_qp];
    case XMOMENTUM:
        return _alpha*_area[_qp]*(density*_vel*_vel+_pressure)*_normals[_qp](0) * _test[_i][_qp];
    case ENERGY:
        return _alpha*_area[_qp]*_vel*(tot_energy+_pressure)*_normals[_qp](0) * _test[_i][_qp];
    default:
    mooseError("The equation with name: \"" << _eqn_name << "\" is not supported in the \"EelNeumannBC\" type of boundary condition.");
    }
}

Real
EelNeumannBC::computeQpJacobian()
{
    return 0.;
}

Real
EelNeumannBC::computeQpOffDiagJacobian(unsigned _jvar)
{
    return 0.;
}
