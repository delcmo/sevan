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
This function computes the fluid internal energy 'rhoe' from the conservative variables. It is dimension agnostic.
**/
#include "InternalEnergyAux.h"

template<>
InputParameters validParams<InternalEnergyAux>()
{
  InputParameters params = validParams<AuxKernel>();
    // Coupled variables
    params.addRequiredCoupledVar("alrhoA", "fluid density: alpha*rho*A");
    params.addRequiredCoupledVar("alrhouA_x", "fluid x momentum component");
    params.addCoupledVar("alrhouA_y", "fluid x momentum component");
    params.addCoupledVar("alrhouA_z", "fluid x momentum component");
    params.addRequiredCoupledVar("alrhoEA", "alpha*rho*E*A");
    // Coupled aux variables
    params.addRequiredCoupledVar("vf_liquid","void fraction of the liquid");
    params.addRequiredCoupledVar("area", "area");
    params.addParam<bool>("isLiquid", true,"is the flud liquid");
  return params;
}

InternalEnergyAux::InternalEnergyAux(const std::string & name, InputParameters parameters) :
    AuxKernel(name, parameters),
    // Coupled variables:
    _alrhoA(coupledValue("alrhoA")),
    _alrhouA_x(coupledValue("alrhouA_x")),
    _alrhouA_y(_mesh.dimension()>=2 ? coupledValue("alrhouA_y") : _zero),
    _alrhouA_z(_mesh.dimension()==3 ? coupledValue("alrhouA_z") : _zero),
    _alrhoEA(coupledValue("alrhoEA")),
    // Aux variables:
    _alpha_liq(coupledValue("vf_liquid")),
    _area(coupledValue("area")),
    // Parameters:
    _isLiquid(getParam<bool>("isLiquid"))
{}

Real
InternalEnergyAux::computeValue()
{
    // Compute the phase void fraction:
    Real _alpha = _isLiquid ? _alpha_liq[_qp] : 1.-_alpha_liq[_qp];
    
    // Compute density, norm of velocity and total energy:
    Real _rhoE = _alrhoEA[_qp] / (_area[_qp]*_alpha);
    Real _vel_x = _alrhouA_x[_qp] / _alrhoA[_qp];
    Real _vel_y = _alrhouA_y[_qp] / _alrhoA[_qp];
    Real _vel_z = _alrhouA_z[_qp] / _alrhoA[_qp];
    Real _rho = _alrhoA[_qp] / (_area[_qp]*_alpha);
    Real _norm_vel2 = _vel_x*_vel_x + _vel_y*_vel_y + _vel_z*_vel_z;
    
    // Return internal energy:
    return _rhoE - 0.5*_rho*_norm_vel2;
}
