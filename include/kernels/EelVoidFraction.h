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

#ifndef EELVOIDFRACTION_H
#define EELVOIDFRACTION_H

#include "Kernel.h"

class EelVoidFraction;

template<>
InputParameters validParams<EelVoidFraction>();
class EelVoidFraction : public Kernel
{
public:

  EelVoidFraction(const std::string & name,
             InputParameters parameters);

protected:
 
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian( unsigned int jvar );

private:
    // Coupled aux variables
    VariableValue & _pressure_l;
    VariableValue & _pressure_g;
    VariableGradient & _grad_alpha_l;
    VariableValue & _area;
    // Coupled material:
    MaterialProperty<RealVectorValue> & _velI;
    MaterialProperty<Real> & _P_rel;
    MaterialProperty<Real> & _rhoI;
    MaterialProperty<Real> & _Aint;
    MaterialProperty<Real> & _Omega_gas;
};

#endif // VOIDFRACTION_H
