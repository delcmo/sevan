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

#ifndef EELMOMENTUM_H
#define EELMOMENTUM_H

#include "Kernel.h"

// Forward Declarations
class EelMomentum;

template<>
InputParameters validParams<EelMomentum>();

class EelMomentum : public Kernel
{
public:

  EelMomentum(const std::string & name,
             InputParameters parameters);

protected:

  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian( unsigned int jvar );

private:
    // Density variable:
    VariableValue & _alrhoA;
   
    // Velocity aux variables:
    VariableValue & _vel_x;
    VariableValue & _vel_y;
    VariableValue & _vel_z;
    VariableValue & _vel_x_2;
    VariableValue & _vel_y_2;
    VariableValue & _vel_z_2;
    
    // Pressure
    VariableValue & _pressure;
    
    // Area and liquid void fraction:
    VariableValue & _area;
    VariableGradient & _grad_area;
    VariableValue &  _alpha_liq;
    VariableGradient & _grad_alpha_liq;
    
    // Phase related parameters:
    int _component;
    bool _isLiquid;
    
    // Material property: interfacial variables.
    MaterialProperty<Real> & _PI;
    MaterialProperty<RealVectorValue> & _velI;
    
    // Material property: realxation parameters.
    MaterialProperty<Real> & _vel_rel;
    
    // Material property: mass transfer
    MaterialProperty<Real> & _Aint;
    MaterialProperty<Real> & _Omega_gas;
    
    // Material property: friction coefficients
    MaterialProperty<Real> & _wall_friction;
    MaterialProperty<Real> & _interf_friction;
    
    // Gravity vector:
    RealVectorValue _gravity;

};

#endif // EELMOMENTUM_H
