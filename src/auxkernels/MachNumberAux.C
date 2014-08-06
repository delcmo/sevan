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
This function compute the Mach number. It is dimension agnostic.
**/
#include "MachNumberAux.h"

template<>
InputParameters validParams<MachNumberAux>()
{
  InputParameters params = validParams<AuxKernel>();
    // Conservative coupled variables:
    params.addRequiredCoupledVar("alrhoA", "rhoA");
    params.addRequiredCoupledVar("alrhouA_x", "x component of momentum");
    params.addCoupledVar("alrhouA_y", "y component of momentum");
    params.addCoupledVar("alrhouA_z", "z component of momentum");
    // Aux coupled variables:
    params.addRequiredCoupledVar("vf_liquid","void fraction of the liquid");
    params.addRequiredCoupledVar("pressure", "pressure");
    params.addRequiredCoupledVar("area", "area");
    params.addRequiredParam<UserObjectName>("eos", "Equation of state");
    // Parameter:
    params.addParam<bool>("isLiquid", true,"is the flud liquid");
  return params;
}

MachNumberAux::MachNumberAux(const std::string & name, InputParameters parameters) :
    AuxKernel(name, parameters),
    // Coupled variables
    _alrhoA(coupledValue("alrhoA")),
    _alrhouA_x(coupledValue("alrhouA_x")),
    _alrhouA_y(_mesh.dimension()>=2 ? coupledValue("alrhouA_y") : _zero),
    _alrhouA_z(_mesh.dimension()==3 ? coupledValue("alrhouA_z") : _zero),
    // Aux coupled variables:
    _alpha_liq(coupledValue("vf_liquid")),
    _pressure(coupledValue("pressure")),
    _area(coupledValue("area")),
    // User Objects for eos
    _eos(getUserObject<EquationOfState>("eos")),
    // Parameters:
    _isLiquid(getParam<bool>("isLiquid"))
{}

Real
MachNumberAux::computeValue()
{
    // Compute the phase void fraction:
    Real _alpha = _isLiquid ? _alpha_liq[_qp] : 1.-_alpha_liq[_qp];
    
    // Compute the norm of velocity and density:
    Real _rho = _alrhoA[_qp] / (_area[_qp]*_alpha);
    Real _vel_x = _alrhouA_x[_qp] / _alrhoA[_qp];
    Real _vel_y = _alrhouA_y[_qp] / _alrhoA[_qp];
    Real _vel_z = _alrhouA_z[_qp] / _alrhoA[_qp];
    Real _norm_vel2 = _vel_x*_vel_x + _vel_y*_vel_y + _vel_z*_vel_z;
    
    // Computes the speed of sounds:
    Real _c2 = _eos.c2_from_p_rho(_rho, _pressure[_qp]);

    // Return the value of the Mach number:
    return std::sqrt(_norm_vel2 / _c2);
}
