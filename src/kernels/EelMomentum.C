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
    params.addCoupledVar("alrhoA", "alpha*rho*A");
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
    params.addParam<RealVectorValue>("gravity", (0., 0., 0.), "gravity vector");
  return params;
}

EelMomentum::EelMomentum(const std::string & name,
                       InputParameters parameters) :
  Kernel(name, parameters),
    // Coupled variable:
    _alrhoA(isCoupled("alrhoA") ? coupledValue("alrhoA") :_zero),
    // Coupled auxilary variables:
    _vel_x(coupledValue("vel_x")),
    _vel_y(_mesh.dimension()>=2 ? coupledValue("vel_y") : _zero),
    _vel_z(_mesh.dimension()==3 ? coupledValue("vel_z") : _zero),
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
    _Aint(getMaterialProperty<Real>("interfacial_area")),
    _Omega_gas(getMaterialProperty<Real>("mass_transfer")),
    // Material property: friction coefficients.
    _wall_friction(_isLiquid ? getMaterialProperty<Real>("wall_friction_liq") : getMaterialProperty<Real>("wall_friction_gas") ),
    _interf_friction(_isLiquid ? getMaterialProperty<Real>("interfacial_friction_liq") : getMaterialProperty<Real>("interfacial_friction_gas") ),
    // Gravity vector:
    _gravity(getParam<RealVectorValue>("gravity"))
{
    if ( _component > 2 )
        mooseError("ERROR: the integer variable 'component' can only take three values: 0, 1 and 2 that correspond to x, y and z momentum components, respectively.");
}

Real EelMomentum::computeQpResidual()
{
    // Sign: the sign of some terms is phase dependent (+ if liquid, - otherwise).
    Real _sign = _isLiquid ? -1. : 1.;

    // Compute void fraction and its derivative of the phase (liquid or vapor):
    Real _alpha = _isLiquid ? _alpha_liq[_qp] : std::fabs(1.-_alpha_liq[_qp]);
    RealVectorValue _grad_alpha =_isLiquid ? _grad_alpha_liq[_qp] : -_grad_alpha_liq[_qp];

    // Velocity vectors: k->phase under consideration, 2->other phase.
    RealVectorValue _vel_k(_vel_x[_qp], _vel_y[_qp], _vel_z[_qp]);
    RealVectorValue _vel_j(_vel_x_2[_qp], _vel_y_2[_qp], _vel_z_2[_qp]);

    // Convection term: _u = alpha*rho*vel*A
    RealVectorValue _convection(_u[_qp]*_vel_x[_qp], _u[_qp]*_vel_y[_qp], _u[_qp]*_vel_z[_qp]);

    // Pressure term: alpha*P*A
    Real _press = _alpha*_pressure[_qp]*_area[_qp];

    // Source terms: alpha*P*dAdx_i and P*A*dalphadx_i
    Real _source_press = _alpha*_pressure[_qp]*_grad_area[_qp](_component);
    Real _source_alpha = _area[_qp]*_PI[_qp]*_grad_alpha(_component);

    // Relaxation term: lambda*A*(u_2 - u_1)
    Real _source_rel = _area[_qp]*_vel_rel[_qp]*(_vel_j(_component)-_vel_k(_component));

    // Mass transfer source term:
    Real _mass = _sign*_area[_qp]*_Aint[_qp]*_velI[_qp](_component)*_Omega_gas[_qp];

    // Wall friction:
    Real alrho = _alrhoA[_qp] / _area[_qp];
    Real wall_frict_term = -_wall_friction[_qp]*alrho*_vel_k.size()*_vel_k(_component)*std::sqrt(libMesh::pi*_area[_qp]);
    
    // Interfacial friction:
    Real rho = alrho / _alpha;
    RealVectorValue relative_vel = _vel_k - _velI[_qp];
    Real interf_frict_term = -0.5*_interf_friction[_qp]*rho*_Aint[_qp]*_area[_qp]*relative_vel.size()*(_vel_k(_component)-_velI[_qp](_component));
    
//    std::cout<<"wall frict="<<_wall_friction[_qp]<<std::endl;
//    std::cout<<"alrho"<<alrho<<std::endl;
//    std::cout<<wall_frict_term<<std::endl;
    
    // Gravity term:
    Real gravity_term = _alrhoA[_qp]*_gravity(_component);

    // Return the kernel value:
    return -(_convection*_grad_test[_i][_qp] + _press*_grad_test[_i][_qp](_component)) - (_source_press+_source_alpha+_source_rel+_mass+wall_frict_term+gravity_term+interf_frict_term)*_test[_i][_qp];
}

Real EelMomentum::computeQpJacobian()
{
  return 0.;
}

Real EelMomentum::computeQpOffDiagJacobian( unsigned int _jvar)
{ 
  return 0.;
}
