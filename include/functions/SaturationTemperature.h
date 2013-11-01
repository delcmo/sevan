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
/* This function compute the saturation temperature from a pressure.*/
#ifndef SATURATIONTEMPERATURE_H
#define SATURATIONTEMPERATURE_H

#include "Function.h"
#include "EquationOfState.h"

class SaturationTemperature;

template<>
InputParameters validParams<SaturationTemperature>();

class SaturationTemperature : public Function
{
public:
  SaturationTemperature(const std::string & name, InputParameters parameters);

  virtual Real value(Real _pressure, const Point & p);

protected:
    // Equation of state:
    const EquationOfState & _eos_liq;
    const EquationOfState & _eos_gas;
};

#endif //SATURATIONTEMPERATURE_H
