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
/**
This function computes the pressure. It is dimension agnostic.
**/
#include "PressureAux.h"

template<>
InputParameters validParams<PressureAux>()
{
  InputParameters params = validParams<AuxKernel>();
    // Conservative variables:
    params.addRequiredCoupledVar("alrhoA", "alpha*rho*A");
    params.addRequiredCoupledVar("alrhouA_x", "alpha*rho*u*A_x");
    params.addCoupledVar("alrhouA_y", "alpha*rho*u*A_y");
    params.addCoupledVar("alrhouA_z", "alpha*rho*u*A_z");
    params.addRequiredCoupledVar("alrhoEA", "alpha*rho*E*A");
    // Aux variables:
    params.addRequiredCoupledVar("vf_liquid", "liquid void fraction");
    params.addRequiredCoupledVar("area", "area");
    // Userobject:
    params.addRequiredParam<UserObjectName>("eos", "Equation of state");
    // Parameters:
    params.addParam<bool>("isLiquid", true, "is the fluid liquid or not");
  return params;
}

PressureAux::PressureAux(const std::string & name, InputParameters parameters) :
    AuxKernel(name, parameters),
  // Coupled variables
    _alrhoA(coupledValue("alrhoA")),
    _alrhouA_x(coupledValue("alrhouA_x")),
    _alrhouA_y(_mesh.dimension()>=2 ? coupledValue("alrhouA_y"): _zero),
    _alrhouA_z(_mesh.dimension()==3 ? coupledValue("alrhouA_z"): _zero),
    _alrhoEA(coupledValue("alrhoEA")),
  // Aux variables:
    _alpha_liq(coupledValue("vf_liquid")),
    _area(coupledValue("area")),
  // User Objects for eos
    _eos(getUserObject<EquationOfState>("eos")),
    // Parameter:
    _isLiquid(getParam<bool>("isLiquid"))

{}

Real
PressureAux::computeValue()
{
    // Compute the phase void fraction:
    Real _alpha = _isLiquid ? _alpha_liq[_qp] : 1.-_alpha_liq[_qp];
    
    // Computes the density, the norm of the velocity and the total energy:
    Real _rho = _alrhoA[_qp] / (_area[_qp]*_alpha);
    Real _rhoE = _alrhoEA[_qp] / (_area[_qp]*_alpha);
    Real _vel_x = _alrhouA_x[_qp] / _alrhoA[_qp];
    Real _vel_y = _alrhouA_y[_qp] / _alrhoA[_qp];
    Real _vel_z = _alrhouA_z[_qp] / _alrhoA[_qp];
    Real _norm_vel = std::sqrt(_vel_x*_vel_x + _vel_y*_vel_y + _vel_z*_vel_z);
    
    // Computes the pressure
    //std::cout<<"pressure="<<_eos.pressure(_rho, _norm_vel, _rhoE)<<std::endl;
    return _eos.pressure(_rho, _norm_vel, _rhoE);
}
