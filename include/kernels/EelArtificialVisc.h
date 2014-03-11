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

#ifndef EELARTIFICIALVISC_H
#define EELARTIFICIALVISC_H

#include "Kernel.h"
#include "EquationOfState.h"

// Forward Declarations
class EelArtificialVisc;

template<>
InputParameters validParams<EelArtificialVisc>();

class EelArtificialVisc : public Kernel
{
public:

  EelArtificialVisc(const std::string & name,
             InputParameters parameters);

protected:

  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian(unsigned int _jvar);
    
private:
    // Equations types
    enum EquationType
    {
        VOID_FRACTION = 0,
        CONTINUITY = 1,
        XMOMENTUM = 2,
        YMOMENTUM = 3,
        ZMOMENTUM = 4,
        ENERGY = 5
    };
    // Diffusion types
    enum DiffusionType
    {
        ENTROPY = 0,
        PARABLOIC = 1
    };
    // Diffusion name
    std::string _equ_name;
    std::string _diff_name;
    
    // Diffusion type
    MooseEnum _equ_type;
    MooseEnum _diff_type;
    
    // Boolean for phase
    bool _isLiquid;
    
    // Coupled aux variables:
    VariableValue & _rho;
    VariableValue & _pressure;
    VariableGradient & _grad_rho;
    VariableGradient & _grad_press;
    VariableValue & _vel_x;
    VariableValue & _vel_y;
    VariableValue & _vel_z;
    VariableGradient & _grad_vel_x;
    VariableGradient & _grad_vel_y;
    VariableGradient & _grad_vel_z;
    VariableValue & _rhoe;
    VariableGradient & _grad_rhoe;
    VariableValue & _area;
    VariableValue & _alpha_liq;
    VariableGradient & _grad_alpha_liq;
    
    // Material property: viscosity coefficient.
    MaterialProperty<Real> & _mu_liq;
    MaterialProperty<Real> & _mu_gas;
    MaterialProperty<Real> & _kappa_liq;
    MaterialProperty<Real> & _kappa_gas;
    MaterialProperty<Real> & _beta;
    
    // Equation of state:
    const EquationOfState & _eos;
};

#endif // EELARTIFICIALVISC_H
