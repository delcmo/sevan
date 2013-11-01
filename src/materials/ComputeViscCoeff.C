#include "ComputeViscCoeff.h"

template<>
InputParameters validParams<ComputeViscCoeff>()
{
  InputParameters params = validParams<Material>();
    // Viscosity type:
    params.addParam<std::string>("viscosity_name", "FIRST_ORDER", "Name of the viscosity definition to use: set to FIRST_ORDER by default.");
    // Boolean for phase:
    params.addParam<bool>("isLiquid", true, "the phase is liquid or not.");
    // Aux variables:
    params.addRequiredCoupledVar("velocity_x", "x component of the velocity");
    params.addCoupledVar("velocity_y", "y component of the velocity");
    params.addCoupledVar("velocity_z", "z component of the velocity");
    params.addRequiredCoupledVar("pressure", "pressure of the fluid");
    params.addRequiredCoupledVar("pressure2", "pressure of the other phase");
    params.addRequiredCoupledVar("density", "density of the fluid: rho");
    params.addCoupledVar("jump_grad_press", "jump of pressure gradient");
    params.addCoupledVar("jump_grad_dens", "jump of density gradient");
    params.addCoupledVar("jump_grad_alpha", "jump of alpha gradient");
    params.addRequiredCoupledVar("norm_velocity", "norm of the velocity vector");
    params.addCoupledVar("vf_liquid", "liquid void fraction.");
    // Constant parameter:
    params.addParam<double>("Ce", 1., "Coefficient for viscosity");
    // Userobject:
    params.addRequiredParam<UserObjectName>("eos", "Equation of state");
    // PPS names:
    params.addParam<std::string>("velocity_PPS_name", "none", "name of the pps for velocity");
    params.addParam<std::string>("alpha_PPS_name", "none", "name of the pps for alpha");
    return params;
}

ComputeViscCoeff::ComputeViscCoeff(const std::string & name, InputParameters parameters) :
    Material(name, parameters),
    // Declare viscosity types
    _visc_name(getParam<std::string>("viscosity_name")),
    _visc_type("LAPIDUS, FIRST_ORDER, FIRST_ORDER_MACH, ENTROPY, INVALID", "INVALID"),
    // Boolean for phase:
    _isLiquid(getParam<bool>("isLiquid")),
    // Liquid void fraction:
    _alpha_l(_isLiquid ? coupledValue("vf_liquid") : _zero),
    _alpha_l_old(_isLiquid ? coupledValueOld("vf_liquid") : _zero),
    _grad_alpha_l(_isLiquid ? coupledGradient("vf_liquid") : _grad_zero),
    _grad_alpha_l_old(_isLiquid ? coupledGradientOld("vf_liquid") : _grad_zero),
    // Velocity variables:
    _vel_x(coupledValue("velocity_x")),
    _vel_y(_dim>=2 ? coupledValue("velocity_y") : _zero),
    _vel_z(_dim==3 ? coupledValue("velocity_z") : _zero),
    _vel_x_old(coupledValueOld("velocity_x")),
    _vel_y_old(_dim>=2 ? coupledValueOld("velocity_y") : _zero),
    _vel_z_old(_dim==3 ? coupledValueOld("velocity_z") : _zero),
    _grad_vel_x(coupledGradient("velocity_x")),
    // Pressure:
    _pressure(coupledValue("pressure")),
    _pressure_old(coupledValueOld("pressure")),
    _pressure2(isCoupled("pressure2") ? coupledValue("pressure2") : coupledValue("pressure")),
    _grad_press(coupledGradient("pressure")),
    _grad_press_old(coupledGradientOld("pressure")),
    // Density:
    _rho(coupledValue("density")),
    _rho_old(coupledValueOld("density")),
    _grad_rho(coupledGradient("density")),
    _grad_rho_old(coupledGradientOld("density")),
    // Jump of pressure, density and alpha gradients:
    _jump_grad_press(isCoupled("jump_grad_press") ? coupledValue("jump_grad_press") : _zero),
    _jump_grad_dens(isCoupled("jump_grad_dens") ? coupledValue("jump_grad_dens") : _zero),
    _jump_grad_alpha(isCoupled("jump_grad_alpha") ? coupledValue("jump_grad_alpha") : _zero),
    // Norm of velocity vector:
    _norm_vel(coupledValue("norm_velocity")),
    // Declare material properties used in mass, momentum and energy equations:
    _mu(_isLiquid ? declareProperty<Real>("mu_liq") : declareProperty<Real>("mu_gas")),
    _mu_max(_isLiquid ? declareProperty<Real>("mu_max_liq") : declareProperty<Real>("mu_max_gas")),
    _kappa(_isLiquid ? declareProperty<Real>("kappa_liq") : declareProperty<Real>("kappa_gas")),
    _kappa_max(_isLiquid ? declareProperty<Real>("kappa_max_liq") : declareProperty<Real>("kappa_max_gas")),
    // Declare material property used in void fraction equation:
    _beta(_isLiquid ? declareProperty<Real>("beta") : declareProperty<Real>("none")),
    _beta_max(_isLiquid ? declareProperty<Real>("beta_max") : declareProperty<Real>("none")),
    // Get interfacial area
    _PIbar(getMaterialProperty<Real>("average_interfacial_pressure")),
    _Prel(getMaterialProperty<Real>("pressure_relaxation")),
    // Get interfacial velocity
    _velI(getMaterialProperty<RealVectorValue>("interfacial_velocity")),
    //_velI_old(getMaterialPropertyOld<RealVectorValue>("interfacial_velocity")),
    // Get parameter Ce
    _Ce(getParam<double>("Ce")),
    // UserObject:
    _eos(getUserObject<EquationOfState>("eos")),
    // PPS name:
    _velocity_pps_name(getParam<std::string>("velocity_PPS_name")),
    _alpha_pps_name(getParam<std::string>("alpha_PPS_name"))
{
    _visc_type = _visc_name;
    if (_Ce < 0.)
        mooseError("The coefficient Ce has to be positive and is in general not larger than 2 when using LAPIDUS.");
}

ComputeViscCoeff::~ComputeViscCoeff()
{
}

void
ComputeViscCoeff::initQpStatefulProperties()
{
}

void
ComputeViscCoeff::computeQpProperties()
{
    // Determine h (length used in definition of first and second order derivatives):
    Real _h = _current_elem->hmin();
    Real _eps = std::sqrt(std::numeric_limits<Real>::min());
    
    // Compute the first order viscosity and the mach number:
    Real _c = std::sqrt(_eos.c2_from_p_rho(_rho[_qp], _pressure[_qp]));
    
    _mu_max[_qp] = 0.5*_h*(_norm_vel[_qp] + _c);
    _kappa_max[_qp] = 0.5*_h*(_norm_vel[_qp] + _c);
    _beta_max[_qp] = 0.5*_h*_velI[_qp].size();
    Real _Mach = _norm_vel[_qp]/_c;
    
    // Initialyze some variables used in the switch statement:
    Real _Dalpha = 0.; Real _jump_press = 0.; Real _jump_alpha = 0.; Real _jump_mu = 0.; Real _jump_kappa = 0.;
    Real _D_P_stt = 0.; Real _D_P = 0.;
    Real _D_rho_stt = 0.; Real _D_rho = 0.;
    Real _max_vel_pps = 0.; Real _alpha_pps = 0.;
    Real _kappa_e = 0.; Real _mu_e = 0.; Real _residual = 0.; Real _norm_kappa = 0.; Real _norm_mu = 0.;
    RealVectorValue _vel(_vel_x[_qp], _vel_y[_qp], _vel_z[_qp]);
    RealVectorValue _vel_old(_vel_x_old[_qp], _vel_y_old[_qp], _vel_z_old[_qp]);
    RealVectorValue _l(0., 0., 0.);
    // Switch statement over viscosity type:
    switch (_visc_type) {
        case LAPIDUS:
            if (_t_step == 1) {
                _mu[_qp] = _mu_max[_qp];
                _kappa[_qp] = _kappa_max[_qp];
                _beta[_qp] = _beta_max[_qp];
            }
            else {
                _mu[_qp] = _Ce*_h*_h*std::fabs(_grad_vel_x[_qp](0));
                _kappa[_qp] = _mu[_qp];
                _beta[_qp] = 0.;
            }
            break;
        case FIRST_ORDER:
            _mu[_qp] = _mu_max[_qp];
            _kappa[_qp] = _kappa_max[_qp];
            if (_t_step == 1) {
                _beta[_qp] = _mu_max[_qp];
            }
            else
                _beta[_qp] = _beta_max[_qp];
            break;
        case FIRST_ORDER_MACH:
            _mu[_qp] = _Mach*_Mach*_mu_max[_qp];
            _kappa[_qp] = _kappa_max[_qp];
            _beta[_qp] = _beta_max[_qp];
            break;
        case ENTROPY:
            // Get the pps values:
            _max_vel_pps = std::max(getPostprocessorValueByName(_velocity_pps_name), _eps);
            _alpha_pps = std::max(getPostprocessorValueByName(_alpha_pps_name), _eps);
            //std::cout<<"vel pps="<<_max_vel_pps<<std::endl;
            //std::cout<<"jump="<<_jump_grad_press[_qp]<<std::endl;
            //std::cout<<"Ce="<<_Ce<<std::endl;
            
            // Compute entropy residual for void fraction equations:
            _Dalpha = (_alpha_l[_qp]-_alpha_l_old[_qp])/_dt + _velI[_qp]*0.5*( _grad_alpha_l[_qp] + _grad_alpha_l_old[_qp] );
            _Dalpha = std::fabs(_alpha_l[_qp]*_Dalpha)/_alpha_pps;///_alpha_l[_qp];
            
            // Compute characteristic equation associated to eigenvalue u:
            _D_P_stt = 0.5*(_vel+_vel_old)*0.5*(_grad_press[_qp]+_grad_press_old[_qp]);
            _D_P = (_pressure[_qp]-_pressure[_qp])/_dt + _D_P_stt;
            _D_rho_stt = 0.5*(_vel+_vel_old)*0.5*(_grad_rho[_qp]+_grad_rho_old[_qp]);
            _D_rho =(_rho[_qp]-_rho_old[_qp])/_dt + _D_rho_stt;
            
            // Compute the normalization factor:
            _norm_mu = 0.5*_rho[_qp]*_c*_c;
            _norm_kappa = 0.5*std::min(_rho[_qp]*_c*_max_vel_pps, _rho[_qp]*_c*_c);
            
            // Compute the jump of gradient of pressure:
            _jump_press = _norm_vel[_qp]*_jump_grad_press[_qp];//(_norm_vel[_qp]+_c)*_jump_grad_press[_qp];
            _jump_mu = std::max( _norm_vel[_qp]*_jump_grad_dens[_qp]/_rho[_qp], _jump_press/_norm_mu );
            _jump_kappa = std::max( _norm_vel[_qp]*_jump_grad_dens[_qp]/_rho[_qp], _jump_press/_norm_kappa );
            _jump_alpha = std::fabs(_velI[_qp].size()*_jump_grad_alpha[_qp])/_alpha_pps;//_alpha_l[_qp];
            
            // Entropy viscosity coefficients:
            _residual = std::fabs(_D_P-_c*_c*_D_rho); // + std::fabs(_PIbar[_qp]*_Prel[_qp]*(_pressure[_qp]-_pressure2[_qp]));
            _kappa_e = _Ce*_h*_h*( _residual/_norm_kappa + 4*_jump_kappa );
            _mu_e = _Ce*_h*_h*( _residual/_norm_mu + 4*_jump_mu );
            
            // Compute the viscosity coefficients:
            if (_t_step <= 1) {
                _mu[_qp] = _mu_max[_qp];
                _kappa[_qp] = _kappa_max[_qp];
                _beta[_qp] = _beta_max[_qp];
            }
            else {
                _beta[_qp] = std::min(_beta_max[_qp], _Ce*_h*_h*(_Dalpha+4*_jump_alpha));
                _kappa[_qp] = std::min( _kappa_max[_qp], _kappa_e );
                _mu[_qp] = std::min( _mu_max[_qp], _mu_e );
            }
            break;
        default:
            mooseError("The viscosity type entered in the input file is not implemented.");
            break;
    }
}
