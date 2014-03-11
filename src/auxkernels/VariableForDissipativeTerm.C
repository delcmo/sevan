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
#include "VariableForDissipativeTerm.h"

template<>
InputParameters validParams<VariableForDissipativeTerm>()
{
  InputParameters params = validParams<AuxKernel>();
    // Conservative variables:
    params.addRequiredCoupledVar("alphaA_liq", "liquid void fraction");
    params.addRequiredCoupledVar("alrhoA", "alpha*rho*A");
    params.addRequiredCoupledVar("alrhouA_x", "alpha*rho*u*A_x");
    params.addCoupledVar("alrhouA_y", "alpha*rho*u*A_y");
    params.addCoupledVar("alrhouA_z", "alpha*rho*u*A_z");
    params.addRequiredCoupledVar("alrhoEA", "alpha*rho*E*A");
    // Aux variables:
    params.addRequiredCoupledVar("area", "area");
    // Userobject:
    params.addRequiredParam<UserObjectName>("eos", "Equation of state");
    // Parameters:
    params.addParam<bool>("isLiquid", true, "is the fluid liquid or not");
  return params;
}

VariableForDissipativeTerm::VariableForDissipativeTerm(const std::string & name, InputParameters parameters) :
    AuxKernel(name, parameters),
  // Coupled variables
    _alA_liq(coupledValue("alphaA_liq")),
    _alrhoA(coupledValue("alrhoA")),
    _alrhouA_x(coupledValue("alrhouA_x")),
    _alrhouA_y(_mesh.dimension()>=2 ? coupledValue("alrhouA_y"): _zero),
    _alrhouA_z(_mesh.dimension()==3 ? coupledValue("alrhouA_z"): _zero),
    _alrhoEA(coupledValue("alrhoEA")),
  // Aux variables:
    _area(coupledValue("area")),
  // User Objects for eos
    _eos(getUserObject<EquationOfState>("eos")),
    // Parameter:
    _isLiquid(getParam<bool>("isLiquid"))

{}

Real
VariableForDissipativeTerm::computeValue()
{
    // Compute the phase void fraction:
    Real alpha = (1-(double)_isLiquid)*(1-_alA_liq[_qp]/_area[_qp]) + (double)_isLiquid*_alA_liq[_qp]/_area[_qp];
    
    // Computes the density, the norm of the velocity, the total energy and the internal energy:
    Real rho = _alrhoA[_qp] / (_area[_qp]*alpha);
    Real rhoE = _alrhoEA[_qp] / (_area[_qp]*alpha);
    Real vel_x = _alrhouA_x[_qp] / _alrhoA[_qp];
    Real vel_y = _alrhouA_y[_qp] / _alrhoA[_qp];
    Real vel_z = _alrhouA_z[_qp] / _alrhoA[_qp];
    Real norm_vel = std::sqrt(vel_x*vel_x + vel_y*vel_y + vel_z*vel_z);
    Real rhoe = rhoE - 0.5*rho*norm_vel*norm_vel;
    
    // Computes the pressure
    Real pressure = _eos.pressure(rho, norm_vel, rhoE);
    
    // Return (rho*P) / (P-rho*e):
    return rho*pressure/(pressure+rhoe);
}
