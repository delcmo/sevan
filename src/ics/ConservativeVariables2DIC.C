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

#include "ConservativeVariables2DIC.h"

template<>
InputParameters validParams<ConservativeVariables2DIC>()
{
    InputParameters params = validParams<InitialCondition>();
    params.addRequiredParam<FunctionName>("area", "function to compute the cross section");
    // Initial conditions:
    params.addRequiredParam<Real>("pressure_init_left", "Initial pressure on the left");
    params.addRequiredParam<Real>("pressure_init_right", "Initial pressure on the right");
    params.addRequiredParam<Real>("vel_init_left", "Initial velocity on the left");
    params.addRequiredParam<Real>("vel_init_right", "Inital velocity on the right");
    params.addRequiredParam<Real>("temp_init_left", "Initial value of the temperature");
    params.addRequiredParam<Real>("temp_init_right", "Initial value of the temperature");
    params.addRequiredParam<Real>("alpha_init_left", "Initial value of the LIQUID void fraction");
    params.addRequiredParam<Real>("alpha_init_right", "Initial value of the LIQUID void fraction");
    // Membrane position:
    params.addRequiredParam<Real>("x_point_source", "Position of the point source: x");
    params.addRequiredParam<Real>("y_point_source", "Position of the point source: y");
    params.addParam<Real>("length", 0.05, "distance from the point source.");
    // Equation of state
    params.addRequiredParam<UserObjectName>("eos", "parameters for eos.");
    // Boolean
    params.addParam<bool>("isLiquid", true, "is phase liquid or not?");
    return params;
}

ConservativeVariables2DIC::ConservativeVariables2DIC(const std::string & name,
                     InputParameters parameters) :
    InitialCondition(name, parameters),
    // Function
    _area(getFunction("area")),
	// IC parameters
    _p_left(getParam<Real>("pressure_init_left")),
    _p_right(getParam<Real>("pressure_init_right")),
    _v_left(getParam<Real>("vel_init_left")),
    _v_right(getParam<Real>("vel_init_right")),
    _t_left(getParam<Real>("temp_init_left")),
    _t_right(getParam<Real>("temp_init_right")),
    _alpha_left(getParam<Real>("alpha_init_left")),
    _alpha_right(getParam<Real>("alpha_init_right")),
    // Position of the membrane:
    _x_pt_source(getParam<Real>("x_point_source")),
    _y_pt_source(getParam<Real>("y_point_source")),
    _length(getParam<Real>("length")),
  	// User Objects
    _eos(getUserObject<EquationOfState>("eos")),
    // Boolean:
    _isLiquid(getParam<bool>("isLiquid"))
{}

Real
ConservativeVariables2DIC::value(const Point & p)
{
// Define and compute parameters used to smooth the initial condition if wished
/*Real _a_p, _b_p, _a_vel, _b_vel, _a_t, _b_t;
_a_p = ( _p_left - _p_right) / ( _x1 - _x2 );
_b_p = ( _x1*_p_right - _x2*_p_left ) / ( _x1 - _x2 );
_a_vel = ( _v_left - _v_right) / ( _x1 - _x2 );
_b_vel = ( _x1*_v_right - _x2*_v_left ) / ( _x1 - _x2 );
_a_t = ( _t_left - _t_right) / ( _x1 - _x2 );
_b_t = ( _x1*_t_right - _x2*_t_left ) / ( _x1 - _x2 );*/
// Get the name of the variable this object acts on
std::string _name_var = _var.name();
// Compute the pressure, velocity and temperature values
Real _pressure = 0.;
Real _temp = 0.;
Real _vel = 0.;
Real _alpha = 0.;
  if ( p(0)>=(_x_pt_source-_length) && p(0)<=(_x_pt_source+_length) )
	{
        if ( p(1)>=(_y_pt_source-_length) && p(1)<=(_y_pt_source+_length) ) {
            _pressure = _p_left;
            _vel = _v_left;
            _temp = _t_left;
            _alpha = (1-(double)_isLiquid)*(1-_alpha_left) + (double)_isLiquid*_alpha_left;
        }
        else {
            _pressure = _p_right;
            _vel = _v_right;
            _temp = _t_right;
             _alpha = (1-(double)_isLiquid)*(1-_alpha_right) + (double)_isLiquid*_alpha_right;
        }
	}
// Compute the conservative variables
Real _density = (_pressure + _eos.Pinf()) / (_eos.Cv()*(_eos.gamma()-1)*_temp);
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
