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
This function computes the void fraction from the variables 'alhpaA' and 'A'. It is dimension agnostic.
**/
#include "VoidFractionAux.h"

template<>
InputParameters validParams<VoidFractionAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("alA", "alpha*A");
  params.addRequiredCoupledVar("area", "area");
  return params;
}

VoidFractionAux::VoidFractionAux(const std::string & name, InputParameters parameters) :
    AuxKernel(name, parameters),
    _alA(coupledValue("alA")),
    _area(coupledValue("area"))
{}

Real
VoidFractionAux::computeValue()
{
    Real _alpha = _alA[_qp] / _area[_qp];
    //std::cout<<"alpha="<<_alpha<<std::endl;
    if ( _alpha < 0 )
        return 0.;
    else if ( _alpha > 1 )
        return 1.;
    else
        return _alA[_qp] / _area[_qp];
}
