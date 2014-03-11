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

#ifndef EELMASS_H
#define EELMAS_H

#include "Kernel.h"

class EelMass;

template<>
InputParameters validParams<EelMass>();
class EelMass : public Kernel
{
public:

  EelMass(const std::string & name,
             InputParameters parameters);

protected:
 
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian( unsigned int jvar );

private:
    // Boolean for phase:
    bool _isLiquid;
    // Coupled variables
    VariableValue & _alrhouA_x;
    VariableValue & _alrhouA_y;
    VariableValue & _alrhouA_z;
    // Coupled aux variable:
    VariableValue & _area;
    // Material property:
    MaterialProperty<Real> & _Aint;
    MaterialProperty<Real> & _Omega_gas;
};

#endif // EelMass_H
