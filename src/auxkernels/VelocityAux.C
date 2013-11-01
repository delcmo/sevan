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
This function computes the velocity components from the density 'rhoA' and the momentum components 'rhouA_{x,y or z}'. It is dimension agnostic.
**/
#include "VelocityAux.h"

template<>
InputParameters validParams<VelocityAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("alrhoA", "alpha*rho*A");
  params.addRequiredCoupledVar("alrhouA", "alpha*rho*u*A_{x,y or z}");
  return params;
}

VelocityAux::VelocityAux(const std::string & name, InputParameters parameters) :
    AuxKernel(name, parameters),
    _alrhoA(coupledValue("alrhoA")),
    _alrhouA(coupledValue("alrhouA"))
{}

Real
VelocityAux::computeValue()
{
  return _alrhouA[_qp] / _alrhoA[_qp];
}
