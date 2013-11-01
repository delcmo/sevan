#ifndef COMPUTEVISCCOEFF_H
#define COMPUTEVISCCOEFF_H

#include "Material.h"
#include "MaterialProperty.h"
#include "EquationOfState.h"

//Forward Declarations
class ComputeViscCoeff;

template<>
InputParameters validParams<ComputeViscCoeff>();

class ComputeViscCoeff : public Material
{
public:
  ComputeViscCoeff(const std::string & name, InputParameters parameters);
  virtual ~ComputeViscCoeff();

protected:
  virtual void initQpStatefulProperties();  
  virtual void computeQpProperties();

private:
    // Viscosity types
    enum ViscosityType
    {
        LAPIDUS = 0,
        FIRST_ORDER = 1,
        FIRST_ORDER_MACH = 2,
        ENTROPY = 3
    };
    // Artificial diffusion name
    std::string _visc_name;
    // Aritifical diffusion type
    MooseEnum _visc_type;
    // Boolean for phase
    bool _isLiquid;
    // Liquid void fraction:
    VariableValue & _alpha_l;
    VariableValue & _alpha_l_old;
    VariableGradient & _grad_alpha_l;
    VariableGradient & _grad_alpha_l_old;
    // Coupled aux variables
    VariableValue & _vel_x;
    VariableValue & _vel_y;
    VariableValue & _vel_z;
    VariableValue & _vel_x_old;
    VariableValue & _vel_y_old;
    VariableValue & _vel_z_old;
    VariableGradient & _grad_vel_x;
    // Pressure:
    VariableValue & _pressure;
    VariableValue & _pressure_old;
    VariableValue & _pressure2; // other pressure phase
    VariableGradient & _grad_press;
    VariableGradient & _grad_press_old;
    // Density:
    VariableValue & _rho;
    VariableValue & _rho_old;
    VariableGradient & _grad_rho;
    VariableGradient & _grad_rho_old;
    // Variables for jump:
    VariableValue & _jump_grad_press;
    VariableValue & _jump_grad_dens;
    VariableValue & _jump_grad_alpha;
    // Norm of the velocity:
    VariableValue & _norm_vel;
    // Material properties: viscosity coefficients.
    MaterialProperty<Real> & _mu;
    MaterialProperty<Real> & _mu_max;
    MaterialProperty<Real> & _kappa;
    MaterialProperty<Real> & _kappa_max;
    MaterialProperty<Real> & _beta;
    MaterialProperty<Real> & _beta_max;
    // Interfactial variable:
    MaterialProperty<Real> & _PIbar;
    MaterialProperty<Real> & _Prel;
    // Material property: interfacial velocity.
    MaterialProperty<RealVectorValue> & _velI;
    //MaterialProperty<RealVectorValue> & _velI_old;
    // Multiplicative coefficient for viscosity:
    double _Ce;
    // UserObject: equation of state
    const EquationOfState & _eos;
    // Name of the posprocessors for pressure, velocity and void fraction:
    std::string _velocity_pps_name;
    std::string _alpha_pps_name;
};

#endif //ComputeViscCoeff_H
