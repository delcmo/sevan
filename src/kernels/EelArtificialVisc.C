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

#include "EelArtificialVisc.h"
/**
This function computes the dissipative terms for all of the equations. It is dimension agnostic.
 */
template<>
InputParameters validParams<EelArtificialVisc>()
{
  InputParameters params = validParams<Kernel>();
    // Equation and diffusion names:
    params.addParam<std::string>("equation_name", "INVALID", "Name of the equation.");
    params.addParam<std::string>("diffusion_name", "PARABOLIC", "Name of the diffusion.");
    // Boolean
    params.addParam<bool>("isLiquid", true, "boolean to determine if liquid phase or not");
    // Boolean for HEM
    params.addParam<bool>("isHEM", false, "is HEM model?");
    // Coupled aux variables
    params.addRequiredCoupledVar("density", "density of the fluid");
    params.addRequiredCoupledVar("pressure", "pressure of the fluid");
    params.addRequiredCoupledVar("velocity_x", "x component of the velocity");
    params.addCoupledVar("velocity_y", "y component of the velocity");
    params.addCoupledVar("velocity_z", "z component of the velocity");
    params.addRequiredCoupledVar("internal_energy", "internal energy of the fluid");
    params.addRequiredCoupledVar("area", "area of the geometry");
    params.addRequiredCoupledVar("vf_liquid", "liquid void fraction");
    params.addRequiredParam<UserObjectName>("eos", "Equation of state");
  return params;
}

EelArtificialVisc::EelArtificialVisc(const std::string & name,
                       InputParameters parameters) :
  Kernel(name, parameters),
    // Declare equation types
    _equ_name(getParam<std::string>("equation_name")),
    _diff_name(getParam<std::string>("diffusion_name")),
    _equ_type("VOID_FRACTION, CONTINUITY, XMOMENTUM, YMOMENTUM, ZMOMENTUM, ENERGY, INVALID", "INVALID"),
    _diff_type("ENTROPY, PARABOLIC, INVALID", "INVALID"),
    // Boolean
    _isLiquid(getParam<bool>("isLiquid")),
    // Boolean for HEM:
    _isHEM(getParam<bool>("isHEM")),
    // Coupled auxilary variables
    _rho(coupledValue("density")),
    _pressure(coupledValue("pressure")),
    _grad_rho(coupledGradient("density")),
    _vel_x(coupledValue("velocity_x")),
    _vel_y(_dim>=2 ? coupledValue("velocity_y") : _zero),
    _vel_z(_dim==3 ? coupledValue("velocity_z") : _zero),
    _grad_vel_x(coupledGradient("velocity_x")),
    _grad_vel_y(_dim>=2 ? coupledGradient("velocity_y") : _grad_zero),
    _grad_vel_z(_dim==3 ? coupledGradient("velocity_z") : _grad_zero),
    _grad_rhoe(coupledGradient("internal_energy")),
    _area(coupledValue("area")),
    _alpha_liq(coupledValue("vf_liquid")),
    _grad_alpha_liq(coupledGradient("vf_liquid")),
    // Get material property: viscosity coefficient.
    _mu_liq(getMaterialProperty<Real>("mu_liq")),
    _mu_gas(getMaterialProperty<Real>("mu_gas")),
    _kappa_liq(getMaterialProperty<Real>("kappa_liq")),
    _kappa_gas(getMaterialProperty<Real>("kappa_gas")),
    _beta(getMaterialProperty<Real>("beta")),
    // Equation of state:
    _eos(getUserObject<EquationOfState>("eos"))
{
    _equ_type = _equ_name;
    _diff_type = _diff_name;
}

Real EelArtificialVisc::computeQpResidual()
{
    // Viscosity coefficient: mu
    Real _mu = (1-(double)_isLiquid)*_mu_gas[_qp] + (double)_isLiquid*_mu_liq[_qp];
    Real _kappa = (1-(double)_isLiquid)*_kappa_gas[_qp] + (double)_isLiquid*_kappa_liq[_qp];
    if (_isHEM) _mu = std::min(_mu_liq[_qp], _mu_gas[_qp]);
    
    // Determine if cell is on boundary or not:
    Real _isOnbnd = 1.;
    if (_current_elem->node(_i) == 0 || _current_elem->node(_i) == _mesh.nNodes()-1)
        _isOnbnd = 0.;
    
    // If statement on diffusion type:
    if (_diff_type == 1) {
        return _mu * _grad_u[_qp] * _grad_test[_i][_qp];
    }
    else if (_diff_type == 0) {
    // Phase void fraction:
    Real _alpha = (1-(double)_isLiquid)*(1-_alpha_liq[_qp]) + (double)_isLiquid*_alpha_liq[_qp];
    RealVectorValue _grad_alpha = -(1-(double)_isLiquid)*_grad_alpha_liq[_qp] + (double)_isLiquid*_grad_alpha_liq[_qp];
        
    // Symmetric gradient of velocity:
    TensorValue<Real> _grad_vel_tensor(_grad_vel_x[_qp], _grad_vel_y[_qp], _grad_vel_z[_qp]);
    TensorValue<Real> _grad_vel_tensor_sym = ( _grad_vel_tensor + _grad_vel_tensor.transpose() ) * 0.5 * _rho[_qp] * _mu;
    
    // Compute f = kappa * grad(rho):
    RealVectorValue _f(_kappa*_grad_rho[_qp](0), _kappa*_grad_rho[_qp](1), _kappa*_grad_rho[_qp](2));
        
    // Compute velocity vector and its norm:
    RealVectorValue _vel_vector(_vel_x[_qp], _vel_y[_qp], _vel_z[_qp]);
    Real _norm_vel2 = _vel_vector.size_sq();
    
    // Compute h = kappa * grad(rho*e):
    RealVectorValue _h(_kappa*_grad_rhoe[_qp](0), _kappa*_grad_rhoe[_qp](1), _kappa*_grad_rhoe[_qp](2));
        
    // Compute term for void fraction equation:
    Real _temp = _eos.temperature_from_p_rho(_pressure[_qp], _rho[_qp]);
    Real _s_rho = -_pressure[_qp]/(_rho[_qp]*_rho[_qp]*_temp);
        
    // return the dissipative terms:
    RealVectorValue _row_f_cross_vel;
    RealVectorValue _row_grad_vel_tensor;
            switch (_equ_type) {
                case VOID_FRACTION:
                    //return _isOnbnd * _area[_qp] * _beta[_qp] * _grad_alpha * _grad_test[_i][_qp];
                    return _isOnbnd*_area[_qp]*_alpha*_beta[_qp]*_s_rho*_grad_rho[_qp]*_grad_test[_i][_qp];
                    break;
                case CONTINUITY: // div(kappa grad(rho))
                    return _isOnbnd * _alpha * _area[_qp] * _f * _grad_test[_i][_qp];
                    break;
                case XMOMENTUM:
                    return _isOnbnd * _alpha * _area[_qp] * ( _vel_x[_qp]*_f + _mu*_rho[_qp]*_grad_vel_x[_qp] ) * _grad_test[_i][_qp];
                    break;
                case YMOMENTUM:
                    return _isOnbnd * _alpha * _area[_qp] * ( _vel_y[_qp]*_f + _mu*_rho[_qp]*_grad_vel_y[_qp] ) * _grad_test[_i][_qp];
                    break;
                case ZMOMENTUM:
                    return _isOnbnd * _alpha * _area[_qp] * ( _vel_z[_qp]*_f + _mu*_rho[_qp]*_grad_vel_z[_qp] ) * _grad_test[_i][_qp];
                    break;
                case ENERGY:
                    return _isOnbnd * _alpha * _area[_qp] * ( _h + 0.5*_f*_norm_vel2 + _grad_vel_tensor_sym*_vel_vector )*_grad_test[_i][_qp];
                    break;
                default:
                    mooseError("INVALID equation name.");
            }
    }
    else {
        mooseError("INVALID diffusion type.");
        return 0.;
    }
}

Real EelArtificialVisc::computeQpJacobian()
{
  return 0.;
}

Real EelArtificialVisc::computeQpOffDiagJacobian( unsigned int _jvar)
{
  return 0.*_jvar;
}
