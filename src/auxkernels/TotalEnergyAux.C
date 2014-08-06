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
This function computes the total energy
**/
#include "TotalEnergyAux.h"

template<>
InputParameters validParams<TotalEnergyAux>()
{
  InputParameters params = validParams<AuxKernel>();
    // Conservative variables:
    params.addRequiredCoupledVar("alrhoEA", "alpha*rho*E*A");
    // Aux variable:
    params.addRequiredCoupledVar("vf_liquid","void fraction of the liquid");
    params.addRequiredCoupledVar("area", "area");
    // Parameter:
    params.addParam<bool>("isLiquid", true,"is the flud liquid");
  return params;
}

TotalEnergyAux::TotalEnergyAux(const std::string & name, InputParameters parameters) :
    AuxKernel(name, parameters),
  // Coupled variables
    _alrhoEA( coupledValue("alrhoEA")),
  // Coupled Aux Variables
    _alpha_liq(coupledValue("vf_liquid")),
    _area(coupledValue("area")),
  // Parameters:
    _isLiquid(getParam<bool>("isLiquid"))
{}

Real
TotalEnergyAux::computeValue()
{
    // Compute the phase void fraction:
    Real _alpha = _isLiquid ? _alpha_liq[_qp] : 1.-_alpha_liq[_qp];
    
    // Return the value of the total energy:
  return _alrhoEA[_qp] / (_alpha*_area[_qp]) ;
}
