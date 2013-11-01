#include "EquationOfState.h"
#include "MooseError.h"

template<>
InputParameters validParams<EquationOfState>()
{
  InputParameters params = validParams<UserObject>();
    params.addParam<Real>("gamma", 0, "  gamma");
    params.addParam<Real>("Pinf", 0, "  P infinity");
    params.addParam<Real>("q", 0, "  q coefficient");
    params.addParam<Real>("q_prime", 1., "  q' coefficient");
    params.addParam<Real>("Cv", 0, "heat capacity at constant volume");
    params.addParam<Real>("k", 1., "material conductivity");
    params.addParam<Real>("visc", 1., "phase material viscosity");
    return params;
}

EquationOfState::EquationOfState(const std::string & name, InputParameters parameters) :
  GeneralUserObject(name, parameters),
    _gamma(getParam<Real>("gamma")),
    _Pinf(getParam<Real>("Pinf")),
    _qcoeff(getParam<Real>("q")),
    _qcoeff_prime(getParam<Real>("q_prime")),
    _Cv(getParam<Real>("Cv")),
    _k(getParam<Real>("k")),
    _visc(getParam<Real>("visc"))
{}

EquationOfState::~EquationOfState()
{
  // Destructor, empty
}

void
EquationOfState::destroy()
{
}

Real EquationOfState::pressure(Real rho, Real vel, Real rhoE) const
{
  Real _e = ( rhoE - 0.5 * rho*vel*vel ) / rho;
  return ( (_gamma-1) * ( _e - _qcoeff) * rho - _gamma * _Pinf );
}

Real EquationOfState::rho_from_p_T(Real pressure, Real temperature) const
{
    return ( (pressure + _Pinf) / ((_gamma-1)*_Cv*temperature) );
}

Real EquationOfState::e_from_p_rho(Real pressure, Real rho) const
{
    return ( (pressure + _gamma*_Pinf) / ((_gamma-1)*rho) + _qcoeff );
}

Real EquationOfState::temperature_from_p_rho(Real pressure, Real rho) const
{
    return ( (pressure + _Pinf) / ((_gamma-1)*_Cv*rho) );
}

Real EquationOfState::c2_from_p_rho(Real rho, Real pressure) const
{
    return ( _gamma * ( pressure + _Pinf ) / rho );
}
