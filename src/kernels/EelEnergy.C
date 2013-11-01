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

#include "EelEnergy.h"
/**
This function computes the convective part of the total energy equation.
 */
template<>
InputParameters validParams<EelEnergy>()
{
  InputParameters params = validParams<Kernel>();
    // Conservative variables:
    params.addRequiredCoupledVar("alrhoA", "density");
    params.addRequiredCoupledVar("alrhouA_x", "x component of rhouA");
    params.addCoupledVar("alrhouA_y", "y component of rhouA");
    params.addCoupledVar("alrhouA_z", "z component of rhouA");
    // Aux variables:
    params.addCoupledVar("vel_x_2", "x component of velocity (other phase)");
    params.addCoupledVar("vel_y_2", "y component of velocity (other phase)");
    params.addCoupledVar("vel_z_2", "z component of velocity (other phase)");
    params.addRequiredCoupledVar("pressure_liq", "pressure_liq");
    params.addRequiredCoupledVar("pressure_gas", "pressure_gas");
    params.addRequiredCoupledVar("area", "area");
    params.addRequiredCoupledVar("vf_liquid","liquid void fraction");
    // Parameters:
    params.addParam<bool>("isLiquid", true, "boolean to determine if liquid phase or not");
    // Equation of state:
    params.addRequiredParam<UserObjectName>("eos", "Equation of state");
    return params;
}

EelEnergy::EelEnergy(const std::string & name,
                       InputParameters parameters) :
  Kernel(name, parameters),
    // Boolean
    _isLiquid(getParam<bool>("isLiquid")),
    // Coupled variables:
    _alrhoA(coupledValue("alrhoA")),
    _alrhouA_x(coupledValue("alrhouA_x")),
    _alrhouA_y(_dim>=2 ? coupledValue("alrhouA_y") : _zero),
    _alrhouA_z(_dim==3 ? coupledValue("alrhouA_z") : _zero),
    // Velocity:
    _vel_x_2(isCoupled("vel_x_2") ? coupledValue("vel_x_2") : _zero),
    _vel_y_2(isCoupled("vel_y_2") ? coupledValue("vel_y_2") : _zero),
    _vel_z_2(isCoupled("vel_z_2") ? coupledValue("vel_z_2") : _zero),
    // Pressure: both phases.
    _pressure_l(coupledValue("pressure_liq")),
    _pressure_g(coupledValue("pressure_gas")),
    // Area and liquid void fraction:
    _area(coupledValue("area")),
    _alpha_liq(coupledValue("vf_liquid")),
    _grad_alpha_liq(coupledGradient("vf_liquid")),
    // Equation of state:
    _eos(getUserObject<EquationOfState>("eos")),
    // Material: interfacial variables.
    _Aint(getMaterialProperty<Real>("interfacial_area")),
    _PI(getMaterialProperty<Real>("interfacial_pressure")),
    _PI_bar(getMaterialProperty<Real>("average_interfacial_pressure")),
    _velI(getMaterialProperty<RealVectorValue>("interfacial_velocity")),
    _velI_bar(getMaterialProperty<RealVectorValue>("average_interfacial_velocity")),
    _EI(_isLiquid ? getMaterialProperty<Real>("liquid_interfacial_energy") : getMaterialProperty<Real>("gas_interfacial_energy")),
    _tempI(getMaterialProperty<Real>("interfacial_temperature")),
    // Material: relaxation parameters.
    _P_rel(getMaterialProperty<Real>("pressure_relaxation")),
    _vel_rel(getMaterialProperty<Real>("velocity_relaxation")),
    // Material: mass transfer.
    _Omega(getMaterialProperty<Real>("mass_transfer")),
    // Matearial: heat transfer coefficient:
    _ht(_isLiquid ? getMaterialProperty<Real>("liquid_heat_transfer") : getMaterialProperty<Real>("gas_heat_transfer"))
{
}

Real EelEnergy::computeQpResidual()
{
    // Sign: the sign of some terms is phase dependent (+ if liquid, - otherwise).
    Real _sign = -(1-(double)_isLiquid) + (double)_isLiquid;
    
    // Compute void fraction and its derivative of the phase (liquid or vapor):
    Real _alpha = (1-(double)_isLiquid)*(1-_alpha_liq[_qp]) + (double)_isLiquid*_alpha_liq[_qp];
    RealVectorValue _grad_alpha = -(1-(double)_isLiquid)*_grad_alpha_liq[_qp] + (double)_isLiquid*_grad_alpha_liq[_qp];
    
    // Velocity vectors: 1->phase under consideration, 2->other phase.
    RealVectorValue _vel_1(_alrhouA_x[_qp]/_alrhoA[_qp], _alrhouA_y[_qp]/_alrhoA[_qp], _alrhouA_z[_qp]/_alrhoA[_qp]);
    RealVectorValue _vel_2(_vel_x_2[_qp], _vel_y_2[_qp], _vel_z_2[_qp]);
    
    // Set the pressure:
    Real _pressure = (1-(double)_isLiquid)*_pressure_g[_qp] + (double)_isLiquid*_pressure_l[_qp];
    
    // Compute convective part of the energy equation:
    RealVectorValue _conv;
    _conv(0) = _alrhouA_x[_qp] * ( _u[_qp] + _alpha*_pressure*_area[_qp] ) / _alrhoA[_qp];
    _conv(1) = _alrhouA_y[_qp] * ( _u[_qp] + _alpha*_pressure*_area[_qp] ) / _alrhoA[_qp];
    _conv(2) = _alrhouA_z[_qp] * ( _u[_qp] + _alpha*_pressure*_area[_qp] ) / _alrhoA[_qp];
    
    // Compute void fraction source term:
    Real _source_alpha = _area[_qp]*_PI[_qp]*_velI[_qp]*_grad_alpha_liq[_qp];
    //std::cout<<"alpha="<<_source_alpha<<std::endl;
    // Velocity relaxation source term:
    Real _source_vel_rel = _area[_qp]*_velI_bar[_qp]*_vel_rel[_qp]*(_vel_2 - _vel_1);
    
    // Pressure relaxation source term:
    Real _source_press_rel = _sign*_area[_qp]*_PI_bar[_qp]*_P_rel[_qp]*(_pressure_g[_qp]-_pressure_l[_qp]);
    //std::cout<<"press_rel="<<_source_press_rel<<std::endl;
    // Mass transfer source term:
    Real _mass = _sign*_area[_qp]*_Aint[_qp]*_EI[_qp]*_Omega[_qp];
    //std::cout<<"mass="<<_mass<<std::endl;
    
    // Heat transfer source term:
    Real _rho = _alrhoA[_qp]/(_alpha*_area[_qp]);
    Real _temp_phase = _eos.temperature_from_p_rho(_pressure, _rho);
    Real _source_ht = _area[_qp]*_Aint[_qp]*_ht[_qp]*(_tempI[_qp]-_temp_phase);
    
    // Total source term:
    Real _source = _source_alpha + _source_vel_rel + _source_press_rel + _mass + _source_ht;
    
    /// Returns the residual
    return -_conv * _grad_test[_i][_qp] - _source * _test[_i][_qp];
}

Real EelEnergy::computeQpJacobian()
{
    return 0;
}

Real EelEnergy::computeQpOffDiagJacobian( unsigned int _jvar)
{
    return 0;
}
