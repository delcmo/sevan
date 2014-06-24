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
    params.addParam<RealVectorValue>("gravity", (0., 0., 0.), "gravity vector");
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
    _alrhouA_y(_mesh.dimension()>=2 ? coupledValue("alrhouA_y") : _zero),
    _alrhouA_z(_mesh.dimension()==3 ? coupledValue("alrhouA_z") : _zero),
    // Velocity:
    _vel_x_2(isCoupled("vel_x_2") ? coupledValue("vel_x_2") : _zero),
    _vel_y_2(isCoupled("vel_y_2") ? coupledValue("vel_y_2") : _zero),
    _vel_z_2(isCoupled("vel_z_2") ? coupledValue("vel_z_2") : _zero),
    // Pressure: both phases.
    _pressure_l(coupledValue("pressure_liq")),
    _pressure_g(coupledValue("pressure_gas")),
    // Area and liquid void fraction:
    _area(coupledValue("area")),
    _grad_area(coupledGradient("area")),
    _alpha_liq(coupledValue("vf_liquid")),
    _grad_alpha_liq(coupledGradient("vf_liquid")),
    // Equation of state:
    _eos(getUserObject<EquationOfState>("eos")),
    // Gravity vector:
    _gravity(getParam<RealVectorValue>("gravity")),
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
    _Omega_gas(getMaterialProperty<Real>("mass_transfer")),
    // Matearial: heat transfer coefficient:
    _interf_ht(_isLiquid ? getMaterialProperty<Real>("liquid_heat_transfer") : getMaterialProperty<Real>("gas_heat_transfer")),
    // Matearial: wall heat transfer coefficient:
    _wall_ht(_isLiquid ? getMaterialProperty<Real>("wall_heat_transfer_liq") : getMaterialProperty<Real>("wall_heat_transfer_gas")),
    // Material: wall temperature.
    _wall_temp(getMaterialProperty<Real>("wall_temperature"))
{
}

Real EelEnergy::computeQpResidual()
{
    // Sign: the sign of some terms is phase dependent (+ if liquid, - otherwise).
    Real _sign = _isLiquid ? -1. : 1.;
    
    // Compute void fraction and its derivative of the phase (liquid or vapor):
    Real _alpha = _isLiquid ? _alpha_liq[_qp] : 1.-_alpha_liq[_qp];
    RealVectorValue _grad_alpha =_isLiquid ? _grad_alpha_liq[_qp] : -_grad_alpha_liq[_qp];
    
    // Velocity vectors: 1->phase under consideration, 2->other phase.
    RealVectorValue _vel_k(_alrhouA_x[_qp]/_alrhoA[_qp], _alrhouA_y[_qp]/_alrhoA[_qp], _alrhouA_z[_qp]/_alrhoA[_qp]);
    RealVectorValue _vel_j(_vel_x_2[_qp], _vel_y_2[_qp], _vel_z_2[_qp]);
    
    // Set the pressure:
    Real _pressure = _isLiquid ? _pressure_l[_qp] : _pressure_g[_qp];
    
    // Compute convective part of the energy equation:
    RealVectorValue _conv;
    _conv(0) = _alrhouA_x[_qp] * ( _u[_qp] + _alpha*_pressure*_area[_qp] ) / _alrhoA[_qp];
    _conv(1) = _alrhouA_y[_qp] * ( _u[_qp] + _alpha*_pressure*_area[_qp] ) / _alrhoA[_qp];
    _conv(2) = _alrhouA_z[_qp] * ( _u[_qp] + _alpha*_pressure*_area[_qp] ) / _alrhoA[_qp];
    
    // Compute void fraction source term:
    Real _source_alpha = _area[_qp]*_PI[_qp]*_velI[_qp]*_grad_alpha;
    
    // Velocity relaxation source term:
    Real _source_vel_rel = _area[_qp]*_velI_bar[_qp]*_vel_rel[_qp]*(_vel_j - _vel_k);
    
    // Pressure relaxation source term:
    Real _source_press_rel = _sign*_area[_qp]*_PI_bar[_qp]*_P_rel[_qp]*(_pressure_l[_qp]-_pressure_g[_qp]);
    
    // Mass transfer source term:
    Real _mass = _sign*_area[_qp]*_Aint[_qp]*_EI[_qp]*_Omega_gas[_qp];
    
    // Interfacial heat transfer source term:
    Real rho = _alrhoA[_qp]/(_alpha*_area[_qp]);
    Real temp_phase = _eos.temperature_from_p_rho(_pressure, rho);
    Real source_interf_ht = _area[_qp]*_Aint[_qp]*_interf_ht[_qp]*(_tempI[_qp]-temp_phase);
    
    // Wall heat source term:
    Real wall_area = std::sqrt(4*libMesh::pi*_area[_qp]+_grad_area[_qp].size_sq());
    Real source_wall_ht = _alpha*_wall_ht[_qp]*(_wall_temp[_qp]-temp_phase)*wall_area;
    
    // Gravity work:
    Real gravity_work = _alrhoA[_qp]*_vel_k*_gravity;
    
//    std::cout<<"gravity="<<_gravity<<std::endl;
//    std::cout<<"wall temp="<<_wall_temp[_qp]<<std::endl;
//    std::cout<<"wall ht="<<_wall_ht[_qp]<<std::endl;
    
    // Total source term:
    Real _source = _source_alpha + _source_vel_rel + _source_press_rel + _mass + source_interf_ht + source_wall_ht + gravity_work;
    
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
