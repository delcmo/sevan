/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "ConservativeVariables1DXIC.h"

template<>
InputParameters validParams<ConservativeVariables1DXIC>()
{
  InputParameters params = validParams<InitialCondition>();
    params.addRequiredParam<FunctionName>("area", "function to compute the cross section");
    // Initial conditions:
    params.addRequiredParam<Real>("pressure_init_left", "Initial pressure on the left");
    params.addRequiredParam<Real>("pressure_init_right", "Initial pressure on the right");
    params.addParam<Real>("pressure_init_left_gas", "Initial GAS pressure on the left");
    params.addParam<Real>("pressure_init_right_gas", "Initial GAS pressure on the right");
    params.addRequiredParam<Real>("vel_init_left", "Initial velocity on the left");
    params.addRequiredParam<Real>("vel_init_right", "Inital velocity on the right");
    params.addParam<Real>("vel_init_left_gas", "Initial GAS velocity on the left");
    params.addParam<Real>("vel_init_right_gas", "Inital GAS velocity on the right");
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
    // Membrane position:
    params.addParam<Real>("membrane", 0.5, "The value of the membrane");
    params.addParam<Real>("length", 0.01, "To smooth the IC over a given length");
    // Equation of state
    params.addRequiredParam<UserObjectName>("eos", "parameters for eos.");
    // Boolean
    params.addParam<bool>("isLiquid", true, "is phase liquid or not?");
  return params;
}

ConservativeVariables1DXIC::ConservativeVariables1DXIC(const std::string & name,
                     InputParameters parameters) :
    InitialCondition(name, parameters),
    // Function
    _area(getFunction("area")),
	// IC parameters
    _p_left_liq(getParam<Real>("pressure_init_left")),
    _p_right_liq(getParam<Real>("pressure_init_right")),
    _p_left_gas(isParamValid("pressure_init_left_gas") ? getParam<Real>("pressure_init_left_gas") : getParam<Real>("pressure_init_left")),
    _p_right_gas(isParamValid("pressure_init_right_gas") ? getParam<Real>("pressure_init_right_gas") : getParam<Real>("pressure_init_right")),
    _v_left_liq(getParam<Real>("vel_init_left")),
    _v_right_liq(getParam<Real>("vel_init_right")),
    _v_left_gas(isParamValid("vel_init_left_gas") ? getParam<Real>("vel_init_left_gas") : getParam<Real>("vel_init_left")),
    _v_right_gas(isParamValid("vel_init_right_gas") ? getParam<Real>("vel_init_right_gas") : getParam<Real>("vel_init_right")),
    _t_left_liq(getParam<Real>("temp_init_left")),
    _t_right_liq(getParam<Real>("temp_init_right")),
    _t_left_gas(isParamValid("temp_init_left_gas") ? getParam<Real>("temp_init_left_gas") : getParam<Real>("temp_init_left")),
    _t_right_gas(isParamValid("temp_init_right_gas") ? getParam<Real>("temp_init_right_gas") : getParam<Real>("temp_init_right")),
    _rho_left_liq(getParam<Real>("rho_init_left_liq")),
    _rho_right_liq(getParam<Real>("rho_init_right_liq")),
    _rho_left_gas(getParam<Real>("rho_init_left_gas")),
    _rho_right_gas(getParam<Real>("rho_init_right_gas")),
    _alpha_left(getParam<Real>("alpha_init_left")),
    _alpha_right(getParam<Real>("alpha_init_right")),
    // Position of the membrane:
    _membrane(getParam<Real>("membrane")),
    _length(getParam<Real>("length")),
  	// User Objects
    _eos(getUserObject<EquationOfState>("eos")),
    // Boolean:
    _isLiquid(getParam<bool>("isLiquid"))
{
    if ( isParamValid("temp_init_left") && isParamValid("temp_init_right") )
        _isDensity = false;
    else
        _isDensity = true;
}

Real
ConservativeVariables1DXIC::value(const Point & p)
{
// Define and compute parameters used to smooth the initial condition if wished
Real _x1 = _membrane - 0.5 * _length;
Real _x2 = _x1 + _length;
Real _a_p, _b_p, _a_vel, _b_vel, _a_t, _b_t, _a_rho, _b_rho, _a_al, _b_al;
Real _p_left = _isLiquid ? _p_left_liq : _p_left_gas;
Real _p_right = _isLiquid ? _p_right_liq : _p_right_gas;
_a_p = ( _p_left - _p_right) / ( _x1 - _x2 );
_b_p = ( _x1*_p_right - _x2*_p_left ) / ( _x1 - _x2 );
Real _v_left = _isLiquid ? _v_left_liq : _v_left_gas;
Real _v_right = _isLiquid ? _v_right_liq : _v_right_gas;
_a_vel = ( _v_left - _v_right) / ( _x1 - _x2 );
_b_vel = ( _x1*_v_right - _x2*_v_left ) / ( _x1 - _x2 );
Real _t_left, _t_right, _rho_right, _rho_left;
if (!_isDensity)
{
    _t_left = _isLiquid ? _t_left_liq : _t_left_gas;
    _t_right = _isLiquid ? _t_right_liq : _t_right_gas;
    _a_t = ( _t_left - _t_right) / ( _x1 - _x2 );
    _b_t = ( _x1*_t_right - _x2*_t_left ) / ( _x1 - _x2 );
}
else
{
    _rho_left = _isLiquid ? _rho_left_liq : _rho_left_gas;
    _rho_right = _isLiquid ? _rho_right_liq : _rho_right_gas;
    _a_rho = ( _rho_left - _rho_right) / ( _x1 - _x2 );
    _b_rho = ( _x1*_rho_right - _x2*_rho_left ) / ( _x1 - _x2 );
}
_a_al = ( _alpha_left - _alpha_right) / ( _x1 - _x2 );
_b_al = ( _x1*_alpha_right - _x2*_alpha_left ) / ( _x1 - _x2 );
// Get the name of the variable this object acts on
std::string _name_var = _var.name();
// Compute the pressure, velocity and temperature values
Real _pressure = 0.;
Real _temp = 0.;
Real _rho = 0.;
Real _vel = 0.;
Real _alpha = 0.;
if ( p(0) <= _x1 )
{
    _pressure = _p_left;
    _vel = _v_left;
    if (!_isDensity)
        _temp = _t_left;
    else
        _rho = _rho_left;
    _alpha = _isLiquid ? _alpha_left : 1. - _alpha_left;
}
else if ( p(0) >= _x2 )
{
    _pressure = _p_right;
    _vel = _v_right;
    if (!_isDensity)
        _temp = _t_right;
    else
        _rho = _rho_right;
    _alpha = _isLiquid ? _alpha_right : 1. - _alpha_right;
}
else
{
    _pressure = ( _a_p * p(0) + _b_p );
    _vel = ( _a_vel * p(0) + _b_vel );
    if (!_isDensity)
        _temp = ( _a_t * p(0) + _b_t );
    else
        _rho = ( _a_rho * p(0) + _b_rho );
    _alpha = _isLiquid ? _a_al*p(0)+_b_al : 1-_a_al*p(0)-_b_al;
}
// Compute the conservative variables
Real _density = 0.;
if (!_isDensity)
    _density = (_pressure + _eos.Pinf()) / (_eos.Cv()*(_eos.gamma()-1)*_temp);
else
    _density = _rho;
Real _int_energy = (_pressure+_eos.gamma()*_eos.Pinf())/(_density*(_eos.gamma()-1)) + _eos.qcoeff();
Real _tot_energy = _density*(_int_energy + 0.5*_vel*_vel);
// Value of the area:
    Real _A = _area.value(0., p);
// Return the value of the initial condition. Identify the name of the variable
if (_name_var == "alA_l")
    return ( _alpha*_A);
// Density: rhoA
else if ( _name_var == "alrhoA_l" || _name_var == "alrhoA_g" )
	return ( _alpha*_density*_A );
// Momentum: rhouA and rhovA
else if ( _name_var == "alrhouA_l" || _name_var == "alrhouA_g" )
{
	return ( _alpha*_density*_vel*_A);
}
else if ( _name_var == "alrhovA_l" || _name_var == "alrhovA_g" )
{
	return ( _alpha*_density*_vel*_A);
}
// total energy: rhoEA
else if ( _name_var == "alrhoEA_l" || _name_var == "alrhoEA_g" )
	return ( _alpha*_tot_energy*_A );
else
	return 0;
}
