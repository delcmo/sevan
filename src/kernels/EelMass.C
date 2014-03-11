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

#include "EelMass.h"

/**
This Kernel computes the convection flux of the continuity equation :
rho*u*A where A is the area of the geometry.
*/
template<>
InputParameters validParams<EelMass>()
{
  InputParameters params = validParams<Kernel>();
    params.addRequiredCoupledVar("alrhouA_x", "x component of rhouA");
    params.addCoupledVar("alrhouA_y", "y component of rhouA");
    params.addCoupledVar("alrhouA_z", "z component of rhouA");
    params.addRequiredCoupledVar("area", "area");
    params.addParam<bool>("isLiquid", true, "is liquid phase or not?");
  return params;
}

EelMass::EelMass(const std::string & name,
                       InputParameters parameters) :
  Kernel(name, parameters),
    // Boolean for phase:
    _isLiquid(getParam<bool>("isLiquid")),
    // Coupled aux variables
    _alrhouA_x(coupledValue("alrhouA_x")),
    _alrhouA_y(_mesh.dimension()>=2 ? coupledValue("alrhouA_y") : _zero ),
    _alrhouA_z(_mesh.dimension()==3 ? coupledValue("alrhouA_z") : _zero ),
    // Coupled aux variable:
    _area(coupledValue("area")),
    // Material properties:
    _Aint(getMaterialProperty<Real>("interfacial_area")),
    _Omega_gas(getMaterialProperty<Real>("mass_transfer"))
{}

Real EelMass::computeQpResidual()
{
    // Sign: the sign of some terms is phase dependent (+ if liquid, - otherwise).
    Real _sign = _isLiquid ? -1. : 1.;
    
    // Compute convective part of the continuity equation:
    RealVectorValue _conv(_alrhouA_x[_qp], _alrhouA_y[_qp], _alrhouA_z[_qp]);
    
    // Mass transfer source terms:
    Real _mass = _sign*_area[_qp]*_Omega_gas[_qp]*_Aint[_qp];
    
    // Return the total expression for the continuity equation:
    return -_conv * _grad_test[_i][_qp] - _mass*_test[_i][_qp];
}

Real EelMass::computeQpJacobian()
{
  return ( 0 );
}

Real EelMass::computeQpOffDiagJacobian( unsigned int _jvar)
{ 
    return ( 0 );
}
