#include "TimeStepCFL.h"

template<>
InputParameters validParams<TimeStepCFL>()
{
  InputParameters params = validParams<ElementPostprocessor>();

  // Coupled variables
  params.addRequiredCoupledVar("volume-fraction-liquid", "Liquid volume fraction");
  params.addRequiredCoupledVar("alpha-A-rho-liq", "alpha*rhoA liquid");
  params.addRequiredCoupledVar("alpha-A-rhou-liq", "alpha*rhouA liquid");
  params.addRequiredCoupledVar("alpha-A-rhoE-liq", "alpha*rhouA liquid");
  params.addRequiredCoupledVar("alpha-A-rho-gas", "alpha*rhoA vapor");
  params.addRequiredCoupledVar("alpha-A-rhou-gas", "alpha*rhouA vapor");
  params.addRequiredCoupledVar("alpha-A-rhoE-gas", "alpha*rhouA vapor");
  // Equation of state
  params.addRequiredParam<UserObjectName>("eos-liq", "Liquid equation of state");
  params.addRequiredParam<UserObjectName>("eos-gas", "Gas equation of state");    
  // Parameter
  params.addParam<Real>("cfl", 0.8, "CFL number to supply by the user");

  return params;
}

TimeStepCFL::TimeStepCFL(const std::string & name, InputParameters parameters) :
    ElementPostprocessor(name, parameters),
    // Coupled variables
    _alpha_A_l(coupledValue("volume-fraction-liquid")),
    _alpha_A_rho_l(coupledValue("alpha-A-rho-liq")),
    _alpha_A_rhou_l(coupledValue("alpha-A-rhou-liq")),
    _alpha_A_rhoE_l(coupledValue("alpha-A-rhoE-liq")),
    _alpha_A_rho_g(coupledValue("alpha-A-rho-gas")),
    _alpha_A_rhou_g(coupledValue("alpha-A-rhou-gas")),
    _alpha_A_rhoE_g(coupledValue("alpha-A-rhoE-gas")),
    // Equation of state:
    _eos_l(getUserObject<EquationOfState>("eos-liq")),
    _eos_g(getUserObject<EquationOfState>("eos-gas")),
    // Parameters
    _cfl(getParam<Real>("cfl")),
    _value(0.)
{
}

TimeStepCFL::~TimeStepCFL()
{
}

void
TimeStepCFL::initialize()
{
  _value = std::numeric_limits<Real>::max();
}

void
TimeStepCFL::execute()
{
  // Compute cell size
  Real h_cell = std::pow(_current_elem->volume(), 1./_mesh.dimension());

  // Loop over quadrature points
  for (unsigned qp = 0; qp < _qrule->n_points(); ++qp)
  {
    // Compute local max liquid eigenvalue
    Real rho_l = _alpha_A_rho_l[qp]/_alpha_A_l[qp];
    Real vel_l = _alpha_A_rhou_l[qp]/_alpha_A_rho_l[qp];
    Real rhoE_l = _alpha_A_rhoE_l[qp]/_alpha_A_l[qp];
    Real pressure_l = _eos_l.pressure(rho_l, vel_l, rhoE_l);
    Real eigen_l = std::fabs(vel_l)+std::sqrt(_eos_l.c2_from_p_rho(rho_l, pressure_l));

    // Compute local max gas eigenvalue
    Real alpha_A_g = 1.-_alpha_A_l[qp];
    Real rho_g = _alpha_A_rho_g[qp]/alpha_A_g;
    Real vel_g = _alpha_A_rhou_g[qp]/_alpha_A_rho_g[qp];
    Real rhoE_g = _alpha_A_rhoE_g[qp]/alpha_A_g;
    Real pressure_g = _eos_l.pressure(rho_g, vel_g, rhoE_g);
    Real eigen_g = std::fabs(vel_g)+std::sqrt(_eos_g.c2_from_p_rho(rho_g, pressure_g));

    // Compute the local liquid time step
    Real dt_l = _cfl*h_cell/eigen_l;

    // Compute the local gas time step
    Real dt_g = _cfl*h_cell/eigen_g;

    // Return value
    _value = std::min(_value, std::min(dt_g, dt_l));
  }
}

Real
TimeStepCFL::getValue()
{
  _communicator.min(_value);
  return _value;
}

void
TimeStepCFL::threadJoin(const UserObject & uo)
{
  const TimeStepCFL & pps = dynamic_cast<const TimeStepCFL &>(uo);
  _value = std::min(_value, pps._value);
}