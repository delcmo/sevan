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

#ifndef VOIDFRACTIONAUX_H
#define VOIDFRACTIONAUX_H

#include "AuxKernel.h"


//Forward Declarations
class VoidFractionAux;

template<>
InputParameters validParams<VoidFractionAux>();

/**
 * Coupled auxiliary value
 */
class VoidFractionAux : public AuxKernel
{
public:

  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  VoidFractionAux(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeValue();

  bool _implicit;
    VariableValue & _alA;
    VariableValue & _area;
};

#endif //VOIDFRACTIONAUX_H
