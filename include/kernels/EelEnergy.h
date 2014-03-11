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

#ifndef EELENERGY_H
#define EElENERGY_H

#include "Kernel.h"
#include "EquationOfState.h"

// Forward Declarations
class EelEnergy;

template<>
InputParameters validParams<EelEnergy>();

class EelEnergy : public Kernel
{
public:

  EelEnergy(const std::string & name,
             InputParameters parameters);

protected:

  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian( unsigned int _jvar);

private:
    // Phase related parameters:
    bool _isLiquid;
    
    // Coupled variables
    VariableValue & _alrhoA;
    VariableValue & _alrhouA_x;
    VariableValue & _alrhouA_y;
    VariableValue & _alrhouA_z;
    
    // Velocity:
    VariableValue & _vel_x_2;
    VariableValue & _vel_y_2;
    VariableValue & _vel_z_2;
    
    // Pressure:
    VariableValue & _pressure_l;
    VariableValue & _pressure_g;
    
    // Area and liquid void fraction:
    VariableValue & _area;
    VariableGradient & _grad_area;
    VariableValue & _alpha_liq;
    VariableGradient & _grad_alpha_liq;
    
    // Equation of state:
    const EquationOfState & _eos;
    
    // Gravity vector:
    RealVectorValue _gravity;
    
    // Material: interfacial vairables.
    MaterialProperty<Real> & _Aint;
    MaterialProperty<Real> & _PI;
    MaterialProperty<Real> & _PI_bar;
    MaterialProperty<RealVectorValue> & _velI;
    MaterialProperty<RealVectorValue> & _velI_bar;
    MaterialProperty<Real> & _EI;
    MaterialProperty<Real> & _tempI;
    
    // Material: relaxation parameters.
    MaterialProperty<Real> & _P_rel;
    MaterialProperty<Real> & _vel_rel;
    
    // Material property: mass transfer
    MaterialProperty<Real> & _Omega_gas;
    
    // Material property: heat transfer coefficient:
    MaterialProperty<Real> & _interf_ht;
    
    // Material propoerty: wall heat transfer.
    MaterialProperty<Real> & _wall_ht;
    MaterialProperty<Real> & _wall_temp;
};

#endif // EELENERGY_H
