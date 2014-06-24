/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                               */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "EelVoidFraction.h"

/**
This Kernel is used for the liquid void fraction equation. The convective term is not integrated per part because it is non-conservative.
*/
template<>
InputParameters validParams<EelVoidFraction>()
{
  InputParameters params = validParams<Kernel>();
    params.addRequiredCoupledVar("pressure_liq", "liquid pressure");
    params.addRequiredCoupledVar("pressure_gas", "gas pressure");
    params.addRequiredCoupledVar("vf_liquid", "liquid void fraction");
    params.addRequiredCoupledVar("area", "area");
  return params;
}

EelVoidFraction::EelVoidFraction(const std::string & name,
                       InputParameters parameters) :
  Kernel(name, parameters),
  // Coupled aux variables:
    _pressure_l(coupledValue("pressure_liq")),
    _pressure_g(coupledValue("pressure_gas")),
    _grad_alpha_l(coupledGradient("vf_liquid")),
    _area(coupledValue("area")),
  // Material property:
    _velI(getMaterialProperty<RealVectorValue>("interfacial_velocity")),
    _P_rel(getMaterialProperty<Real>("pressure_relaxation")),
    _rhoI(getMaterialProperty<Real>("interfacial_density")),
    _Aint(getMaterialProperty<Real>("interfacial_area")),
    _Omega_gas(getMaterialProperty<Real>("mass_transfer"))
{}

Real EelVoidFraction::computeQpResidual()
{
    // Compute convective part of the liquid void fraction equation:
    Real _conv = _area[_qp]*(_velI[_qp]*_grad_alpha_l[_qp]);
    
    // Relaxation term:
    Real _rel = _area[_qp]*_P_rel[_qp]*(_pressure_l[_qp]-_pressure_g[_qp]);
    
    // Mass trasnfer term:
    Real _mass = _area[_qp]*_Omega_gas[_qp]*_Aint[_qp]/_rhoI[_qp];
    
    // Return the total expression for the liquid void fraction equation:
    return (_conv - _rel + _mass) * _test[_i][_qp];
}

Real EelVoidFraction::computeQpJacobian()
{
  return ( 0 );
}

Real EelVoidFraction::computeQpOffDiagJacobian( unsigned int _jvar)
{ 
    return ( 0 );
}
