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
    
    // Function Mach number
    enum FctOfMachType
    {
        MACH = 0,
        SQRT_MACH = 1,
        FCT_OF_MACH = 2
    };
    
    // Artificial diffusion name
    std::string _visc_name;
    
    // Aritifical diffusion type
    MooseEnum _visc_type;
    
    std::string _fct_of_mach_name;
    MooseEnum _fct_of_mach_type;
    
    // Boolean for phase
    bool _isLiquid;
    bool _isJumpOn;
    
    // Bool for viscosity coefficient:
    bool _isShock;
    bool _areViscEqual;
    
    // Liquid void fraction:
    VariableValue & _alpha_l;
    VariableValue & _alpha_l_old;
    VariableValue & _alpha_l_older;
    VariableGradient & _grad_alpha_l;
    
    // Coupled aux variables
    VariableValue & _vel_x;
    VariableValue & _vel_y;
    VariableValue & _vel_z;
    VariableGradient & _grad_vel_x;
    
    // Pressure:
    VariableValue & _pressure;
    VariableValue & _pressure_old;
    VariableValue & _pressure_older;
    VariableGradient & _grad_press;
    
    // Density:
    VariableValue & _rho;
    VariableValue & _rho_old;
    VariableValue & _rho_older;
    VariableGradient & _grad_rho;
    
    // Variables for jump:
    VariableValue & _jump_grad_press;
    VariableValue & _jump_grad_dens;
    VariableValue & _jump_grad_alpha;
    
    // Area
    VariableValue & _area;
    VariableGradient & _grad_area;
    
    // Material properties: viscosity coefficients.
    MaterialProperty<Real> & _mu;
    MaterialProperty<Real> & _mu_max;
    MaterialProperty<Real> & _kappa;
    MaterialProperty<Real> & _kappa_max;
    MaterialProperty<Real> & _beta;
    MaterialProperty<Real> & _beta_max;
    
    // Material property: interfacial velocity.
    MaterialProperty<RealVectorValue> & _velI;

    // Multiplicative coefficient for viscosity:
    double _Cmax;
    double _Ce;
    double _Cjump;
    double _Calpha;
    
    // Coefficients for 'sigma' function:
    Real _a_coeff;
    Real _Mthres;
    
    // UserObject: equation of state
    const EquationOfState & _eos;
    
    // Name of the posprocessors for rhov2 and void fraction:
    std::string _rhov2_pps_name;
    std::string _alpha_pps_name;
};

#endif //ComputeViscCoeff_H
