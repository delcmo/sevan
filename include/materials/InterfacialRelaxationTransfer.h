#ifndef INTERFACIALRELAXATIONTRANSFER_H
#define INTERFACIALRELAXATIONTRANSFER_H

#include "Material.h"
#include "MaterialProperty.h"
#include "EquationOfState.h"
#include "SaturationTemperature.h"

//Forward Declarations
class InterfacialRelaxationTransfer;

template<>
InputParameters validParams<InterfacialRelaxationTransfer>();

class InterfacialRelaxationTransfer : public Material
{
public:
  InterfacialRelaxationTransfer(const std::string & name, InputParameters parameters);

protected:
  virtual void computeQpProperties();

private:
    // Variable for definition of the interfacial variables:
    enum InterfacialType
    {
        BERRY = 0,
        AMBROSSO = 1,
        LIANG = 2
    };
    std::string _interf_def_name;
    MooseEnum _interf_def;
    
    Real _xi;
    
    // Boolean for mass and heat transfers:
    bool _isMassOn;
    bool _isHeatOn;
    bool _isWallHeatOn;
    bool _isWallFrictOn;
    
    // Boolean for interfacial relaxation parameters:
    bool _isPressRelOn;
    bool _isVelRelOn;
    
    // Coupled aux variables liquid phase
    VariableValue & _vel_x_l;
    VariableValue & _vel_y_l;
    VariableValue & _vel_z_l;
    VariableValue & _pressure_l;
    VariableValue & _rho_l;
    VariableValue & _alpha_l;
    VariableGradient & _grad_alpha_l;
    
    // Coupled aux variables gas phase
    VariableValue & _vel_x_g;
    VariableValue & _vel_y_g;
    VariableValue & _vel_z_g;
    VariableValue & _pressure_g;
    VariableValue & _rho_g;
    
    // Interfacial variables:
    MaterialProperty<Real> & _Aint;
    MaterialProperty<Real> & _PI;
    MaterialProperty<RealVectorValue> & _velI;
    MaterialProperty<Real> & _velI_norm;
    MaterialProperty<Real> & _PI_bar;
    MaterialProperty<RealVectorValue> & _velI_bar;
    MaterialProperty<Real> & _tempI;
    MaterialProperty<Real> & _rhoI;
    MaterialProperty<Real> & _EI_liq;
    MaterialProperty<Real> & _EI_gas;
    
    // Relaxation parameters:
    MaterialProperty<Real> & _P_rel;
    MaterialProperty<Real> & _vel_rel;
    
    // Heat transfer coefficient for liquid and gas phases:
    MaterialProperty<Real> & _ht_liq;
    MaterialProperty<Real> & _ht_gas;
 
    // Mass transfer:
    MaterialProperty<Real> & _Omega_gas;
    
    // Parameters supplied by the user:
    Real _Aint_max;
    Real _wall_heat_liq_value;
    Real _wall_heat_gas_value;
    Real _wall_frict_liq_value;
    Real _wall_frict_gas_value;
    Real _Twall;
    
    // UserObject: equation of state
    const EquationOfState & _eos_liq;
    const EquationOfState & _eos_gas;
    
    // Wall heat transfer for liquid and gas phase
    MaterialProperty<Real> & _wall_ht_liq;
    MaterialProperty<Real> & _wall_ht_gas;
    
    // Friction parameters for liquid and gas phase
    MaterialProperty<Real> & _wall_frict_liq;
    MaterialProperty<Real> & _wall_frict_gas;
    
    // Friction parameters for liquid and gas phase
    MaterialProperty<Real> & _interf_frict_liq;
    MaterialProperty<Real> & _interf_frict_gas;
    
    // Wall temperature:
    MaterialProperty<Real> & _wall_temp;
};

#endif //INTERFACIALRELAXATIONTRANSFER_H