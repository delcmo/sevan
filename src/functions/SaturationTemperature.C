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

#include "SaturationTemperature.h"
/* This function computes the saturation temperature for a pressure. It is based on what is proposed in the DEM paper by Berry and Saurel. */

template<>
InputParameters validParams<SaturationTemperature>()
{
  InputParameters params = validParams<Function>();
    // Equation of state:
    params.addRequiredParam<UserObjectName>("eos_liq", "Liquid equation of state.");
    params.addRequiredParam<UserObjectName>("eos_gas", "Vapor equation of state.");
  return params;
}

SaturationTemperature::SaturationTemperature(const std::string & name, InputParameters parameters) :
    Function(name, parameters),
    // Euqation of state:
    _eos_liq(getUserObject<EquationOfState>("eos_liq")),
    _eos_gas(getUserObject<EquationOfState>("eos_gas"))
{}

Real
SaturationTemperature::value( Real _pressure , const Point & _p)
{
    // Define some parameters used in the local Newton solve: (DEM paper)
    Real _A = (_eos_liq.Cp() - _eos_gas.Cp() + _eos_gas.qcoeff_prime() - _eos_liq.qcoeff_prime()) / (_eos_gas.Cp() - _eos_gas.Cv());
    Real _B = (_eos_liq.qcoeff() - _eos_gas.qcoeff()) / (_eos_gas.Cp() - _eos_gas.Cv());
    Real _C = (_eos_gas.Cp() - _eos_liq.Cp()) / (_eos_gas.Cp() - _eos_gas.Cv());
    Real _D = (_eos_liq.Cp() - _eos_liq.Cv()) / (_eos_gas.Cp() - _eos_gas.Cv());
    
    // Compute the constant residual _R:
    Real _p_term_liq = std::log(_pressure+_eos_liq.Pinf());
    Real _p_term_gas = std::log(_pressure+_eos_gas.Pinf());
    Real _R = _A + _D * _p_term_liq - _p_term_gas;
    
    // Initialyze some values and pick a guess for the temperature:
    Real _f_norm = 1;
    Real _f = 0.0;
    Real _f_prime = 0.0;
    Real _temp = 400.;
    
    // Newton solve:
    while ( std::fabs(_f_norm) > 1e-3)
	{
        _f = _R + _B / _temp + _C * std::log(_temp);
        _f_prime = _C / _temp - _B / (_temp*_temp);
        _temp = _temp - _f / _f_prime;
        _f_norm = _f / _f_prime;
        
	}
    
    // Return the value:
    return _temp;
}

