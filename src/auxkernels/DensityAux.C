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
/**
This function computes the density of the fluid.
**/
#include "DensityAux.h"

template<>
InputParameters validParams<DensityAux>()
{
  InputParameters params = validParams<AuxKernel>();
    params.addRequiredCoupledVar("alrhoA", "alpha*rho*A");
    params.addRequiredCoupledVar("vf_liquid","liquid void fraction");
    params.addRequiredCoupledVar("area", "area");
    params.addParam<bool>("isLiquid", true,"is the fluid liquid or not");
  return params;
}

DensityAux::DensityAux(const std::string & name, InputParameters parameters) :
    AuxKernel(name, parameters),
    _alrhoA(coupledValue("alrhoA")),
    _alpha_liq(coupledValue("vf_liquid")),
    _area(coupledValue("area")),
    _isLiquid(getParam<bool>("isLiquid"))
{}

Real
DensityAux::computeValue()
{
    // Compute the phase void fraction:
    Real _alpha = _isLiquid ? _alpha_liq[_qp] : 1.-_alpha_liq[_qp];
    
    // Return the value of the density:
  return _alrhoA[_qp] / (_area[_qp]*_alpha);
}
