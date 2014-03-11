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

#ifndef VARIABLEFORDISSIPATIVETERM_H
#define VARIABLEFORDISSIPATIVETERM_H

#include "AuxKernel.h"
#include "EquationOfState.h"

//Forward Declarations
class VariableForDissipativeTerm;

template<>
InputParameters validParams<VariableForDissipativeTerm>();

class VariableForDissipativeTerm : public AuxKernel
{
public:

  VariableForDissipativeTerm(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeValue();

    VariableValue & _alA_liq;
    VariableValue & _alrhoA;
    VariableValue & _alrhouA_x;
    VariableValue & _alrhouA_y;
    VariableValue & _alrhouA_z;
    VariableValue & _alrhoEA;
    VariableValue & _area;
    const EquationOfState & _eos;
    bool _isLiquid;
};

#endif //PRESSUREAUX_H
