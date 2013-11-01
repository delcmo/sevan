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

#include "EelMomentum.h"
/**
This function computes the x, y and z momentum equationS. It is dimension agnostic. 
 */
template<>
InputParameters validParams<EelMomentum>()
{
    InputParameters params = validParams<Kernel>();
    params.addRequiredCoupledVar("vel_x", "x component of velocity");
    params.addCoupledVar("vel_y", "y component of velocity");
    params.addCoupledVar("vel_z", "z component of velocity");
    params.addCoupledVar("vel_x_2", "x component of velocity (other phase)");
    params.addCoupledVar("vel_y_2", "y component of velocity (other phase)");
    params.addCoupledVar("vel_z_2", "z component of velocity (other phase)");
    params.addRequiredCoupledVar("pressure", "pressure");
    params.addRequiredCoupledVar("area", "area");
    params.addRequiredCoupledVar("vf_liquid","liquid void fraction");
    params.addParam<int>("component", 0, "component of the momentum equation to compute (0,1,2)->(x,y,z)");
    params.addParam<bool>("isLiquid", true, "boolean to determine if liquid phase or not");
  return params;
}

EelMomentum::EelMomentum(const std::string & name,
                       InputParameters parameters) :
  Kernel(name, parameters),
    // Coupled auxilary variables:
    _vel_x(coupledValue("vel_x")),
    _vel_y(_dim>=2 ? coupledValue("vel_y") : _zero),
    _vel_z(_dim==3 ? coupledValue("vel_z") : _zero),
    _vel_x_2(isCoupled("vel_x_2") ? coupledValue("vel_x_2") : _zero),
    _vel_y_2(isCoupled("vel_y_2") ? coupledValue("vel_y_2") : _zero),
    _vel_z_2(isCoupled("vel_z_2") ? coupledValue("vel_z_2") : _zero),
    _pressure(coupledValue("pressure")),
    _area(coupledValue("area")),
    _grad_area(coupledGradient("area")),
    _alpha_liq(coupledValue("vf_liquid")),
    _grad_alpha_liq(coupledGradient("vf_liquid")),
    // Parameters:
    _component(getParam<int>("component")),
    _isLiquid(getParam<bool>("isLiquid")),
    // Material property: interfacial variables.
    _PI(getMaterialProperty<Real>("interfacial_pressure")),
    _velI(getMaterialProperty<RealVectorValue>("interfacial_velocity")),
    // Material property: relaxation parameters.
    _vel_rel(getMaterialProperty<Real>("velocity_relaxation")),
    // Material property: mass transfer
    _Omega(getMaterialProperty<Real>("mass_transfer"))
{
    if ( _component > 2 )
        mooseError("ERROR: the integer variable 'component' can only take three values: 0, 1 and 2 that correspond to x, y and z momentum components, respectively.");
}

Real EelMomentum::computeQpResidual()
{
  // Sign: the sign of some terms is phase dependent (+ if liquid, - otherwise).
    Real _sign = -(1-(double)_isLiquid) + (double)_isLiquid;
    
  // Compute void fraction and its derivative of the phase (liquid or vapor):
    Real _alpha = (1-(double)_isLiquid)*(1-_alpha_liq[_qp]) + (double)_isLiquid*_alpha_liq[_qp];
    RealVectorValue _grad_alpha = -(1-(double)_isLiquid)*_grad_alpha_liq[_qp] + (double)_isLiquid*_grad_alpha_liq[_qp];
    
  // Velocity vectors: 1->phase under consideration, 2->other phase.
    RealVectorValue _vel_1(_vel_x[_qp], _vel_y[_qp], _vel_z[_qp]);
    RealVectorValue _vel_2(_vel_x_2[_qp], _vel_y_2[_qp], _vel_z_2[_qp]);
    
  // Convection term: _u = alpha*rho*vel*A
    RealVectorValue _convection(_u[_qp]*_vel_x[_qp], _u[_qp]*_vel_y[_qp], _u[_qp]*_vel_z[_qp]);
    
  // Pressure term: alpha*P*A
    Real _press = _alpha*_pressure[_qp]*_area[_qp];
    
  // Source terms: alpha*P*dAdx_i and P*A*dalphadx_i
    Real _source_press = _alpha*_pressure[_qp]*_grad_area[_qp](_component);
    Real _source_alpha = _area[_qp]*_PI[_qp]*_grad_alpha(_component);
    
  // Relaxation term: lambda*A*(u_2 - u_1)
    Real _source_rel = _area[_qp]*_vel_rel[_qp]*(_vel_2(_component)-_vel_1(_component));
    
  // Mass transfer source term:
    Real _mass = _sign*_area[_qp]*_velI[_qp](_component)*_Omega[_qp];
    
  // Return the kernel value:
    return -( _convection*_grad_test[_i][_qp] + _press*_grad_test[_i][_qp](_component) + (_source_press+_source_alpha+_source_rel-_mass)*_test[_i][_qp] );
}

Real EelMomentum::computeQpJacobian()
{
  return 0.;
}

Real EelMomentum::computeQpOffDiagJacobian( unsigned int _jvar)
{ 
  return 0.;
}
