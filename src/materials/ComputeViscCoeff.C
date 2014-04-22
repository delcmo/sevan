#include "ComputeViscCoeff.h"

template<>
InputParameters validParams<ComputeViscCoeff>()
{
  InputParameters params = validParams<Material>();
    // Viscosity type:
    params.addParam<std::string>("viscosity_name", "FIRST_ORDER", "Name of the viscosity definition to use: set to FIRST_ORDER by default.");
    params.addParam<std::string>("function_of_mach", "MACH", "Name of the Mach function to use.");
    // Boolean for phase:
    params.addParam<bool>("isLiquid", true, "the phase is liquid or not.");
    params.addParam<bool>("useVelPps", true, "use velocity pps");
    params.addParam<bool>("usePressPps", false, "use pressure pps");
    params.addParam<bool>("useAlphaPps", false, "use alpha pps");
    // Aux variables:
    params.addRequiredCoupledVar("velocity_x", "x component of the velocity");
    params.addCoupledVar("velocity_y", "y component of the velocity");
    params.addCoupledVar("velocity_z", "z component of the velocity");
    params.addRequiredCoupledVar("pressure", "pressure of the fluid");
    params.addRequiredCoupledVar("density", "density of the fluid: rho");
    params.addCoupledVar("jump_grad_press", "jump of pressure gradient");
    params.addCoupledVar("jump_grad_dens", "jump of density gradient");
    params.addCoupledVar("jump_grad_alpha", "jump of alpha gradient");
    params.addRequiredCoupledVar("norm_velocity", "norm of the velocity vector");
    params.addCoupledVar("vf_liquid", "liquid void fraction.");
    // Constant parameter:
    params.addParam<double>("Ce", 1., "Coefficient for residual");
    params.addParam<double>("Cjump", 1., "Coefficient for jump");
    params.addParam<double>("Calpha", 1., "Coefficient for alpha");
    // Userobject:
    params.addRequiredParam<UserObjectName>("eos", "Equation of state");
    // PPS names:
    params.addParam<std::string>("velocity_PPS_name", "name of the pps for velocity");
    params.addParam<std::string>("pressure_PPS_name", "name of the pps for pressure");
    params.addParam<std::string>("alpha_PPS_name", "name of the pps for alpha");
    return params;
}

ComputeViscCoeff::ComputeViscCoeff(const std::string & name, InputParameters parameters) :
    Material(name, parameters),
    // Declare viscosity types
    _visc_name(getParam<std::string>("viscosity_name")),
    _visc_type("LAPIDUS, FIRST_ORDER, FIRST_ORDER_MACH, ENTROPY, INVALID", "INVALID"),
    // Function Mach number:
    _fct_of_mach_name(getParam<std::string>("function_of_mach")),
    _fct_of_mach_type("MACH, SQRT_MACH, FCT_OF_MACH, INVALID", _fct_of_mach_name),
    // Boolean for phase:
    _isLiquid(getParam<bool>("isLiquid")),
    _useVelPps(getParam<bool>("useVelPps")),
    _usePressPps(getParam<bool>("usePressPps")),
    _useAlphaPps(getParam<bool>("useAlphaPps")),
    // Liquid void fraction:
    _alpha_l(_isLiquid ? coupledValue("vf_liquid") : _zero),
    _alpha_l_old(_isLiquid ? coupledValueOld("vf_liquid") : _zero),
    _alpha_l_older(_isLiquid ? coupledValueOlder("vf_liquid") : _zero),
    _grad_alpha_l(_isLiquid ? coupledGradient("vf_liquid") : _grad_zero),
    // Velocity variables:
    _vel_x(coupledValue("velocity_x")),
    _vel_y(_mesh.dimension()>=2 ? coupledValue("velocity_y") : _zero),
    _vel_z(_mesh.dimension()==3 ? coupledValue("velocity_z") : _zero),
    _grad_vel_x(coupledGradient("velocity_x")),
    // Pressure:
    _pressure(coupledValue("pressure")),
    _pressure_old(coupledValueOld("pressure")),
    _pressure_older(coupledValueOlder("pressure")),
    _grad_press(coupledGradient("pressure")),
    // Density:
    _rho(coupledValue("density")),
    _rho_old(coupledValueOld("density")),
    _rho_older(coupledValueOlder("density")),
    _grad_rho(coupledGradient("density")),
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
    // Get parameter Ce
    _Ce(getParam<double>("Ce")),
    _Cjump(getParam<double>("Cjump")),
    _Calpha(getParam<double>("Calpha")),
    // UserObject:
    _eos(getUserObject<EquationOfState>("eos")),
    // PPS name:
    _velocity_pps_name(getParam<std::string>("velocity_PPS_name")),
    _pressure_pps_name(getParam<std::string>("pressure_PPS_name")),
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
    
    _mu_max[_qp] = 0.5*_h*_norm_vel[_qp];
    _kappa_max[_qp] = 0.5*_h*(_norm_vel[_qp] + _c);
    _beta_max[_qp] = 0.5*_h*_velI[_qp].size();
    Real _Mach = _norm_vel[_qp]/_c;
    Real fct_of_mach = _Mach;
    switch (_fct_of_mach_type) {
        case MACH:
            fct_of_mach = std::min(_Mach, 1.);
            break;
        case SQRT_MACH:
            fct_of_mach = std::min(std::sqrt(_Mach), 1.);
            break;
        case FCT_OF_MACH:
            fct_of_mach = std::min(_Mach*std::sqrt(4+(1.-_Mach*_Mach)*(1.-_Mach*_Mach)) / (1.+_Mach*_Mach),1.);
            break;
        default:
            mooseError("The function with name: \"" << _fct_of_mach_name << "\" is not supported in the \"ComputeViscCoeff\" type of material.");
    }
    
    Real vel_var = _useVelPps ? std::max(getPostprocessorValueByName(_velocity_pps_name), _eps) : _norm_vel[_qp];
    Real press_var = _usePressPps ? std::max(getPostprocessorValueByName(_pressure_pps_name), _eps) : _rho[_qp]*_c*_c;
    Real alpha_var = _useAlphaPps ? getPostprocessorValueByName(_alpha_pps_name) : _alpha_l[_qp];
    
    // Initialyze some variables used in the switch statement:
    Real _Dalpha = 0.; Real _jump = 0.;
    Real _D_P = 0.; Real _jump_alpha = 0.;
    Real _D_stt = 0.; Real _D_rho = 0.;
//    Real _alpha_pps = 0.;
    Real _kappa_e = 0.; Real _mu_e = 0.;
    Real _residual = 0.; Real _norm_kappa = 0.; Real _norm_mu = 0.; Real _norm_alpha = 0.;
    RealVectorValue _vel(_vel_x[_qp], _vel_y[_qp], _vel_z[_qp]);
    RealVectorValue _l(0., 0., 0.);
    Real _weight0 = 0.; Real _weight1 = 0.; Real _weight2 = 0.;
    
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
            // Compute the weigth for BDF2
            _weight0 = (2.*_dt+_dt_old)/(_dt*(_dt+_dt_old));
            _weight1 = -(_dt+_dt_old)/(_dt*_dt_old);
            _weight2 = _dt/(_dt_old*(_dt+_dt_old));
            
            // Compute entropy residual for void fraction equations:
            _D_stt = _velI[_qp]*_grad_alpha_l[_qp];
            _Dalpha = (_weight0*_alpha_l[_qp]+_weight1*_alpha_l_old[_qp]+_weight2*_alpha_l_older[_qp]) + _D_stt;
            
            // Compute the characteristic equation u:
            _D_stt = _vel*_grad_press[_qp];
            _D_P = (_weight0*_pressure[_qp]+_weight1*_pressure_old[_qp]+_weight2*_pressure_older[_qp]) + _D_stt;
            _D_stt = _vel*_grad_rho[_qp];
            _D_rho = (_weight0*_rho[_qp]+_weight1*_rho_old[_qp]+_weight2*_rho_older[_qp]) + _D_stt;
            _residual = std::fabs(_D_P-_c*_c*_D_rho);
            
            // Compute the normalization factor:
            _norm_mu = 0.5*((1.-fct_of_mach)*press_var+fct_of_mach*std::min(vel_var*vel_var, _c*_c));
            _norm_kappa = _norm_mu;
            _norm_alpha = alpha_var;// _alpha_l[_qp];
            
            // Compute the jump of gradient of pressure:
            _jump = _Cjump*_norm_vel[_qp]*std::max( _c*_c*_jump_grad_dens[_qp], _jump_grad_press[_qp]);
            _jump_alpha = _Calpha*std::fabs(_velI[_qp].size()*_jump_grad_alpha[_qp]);
            
            // Entropy viscosity coefficients:
            _kappa_e = _Ce*_h*_h*( _residual + _jump )/_norm_kappa;
            _mu_e = _Ce*_h*_h*( _residual + _jump )/_norm_mu;
            
            // Compute the viscosity coefficients:
            if (_t_step <= 1) {
                _mu[_qp] = _mu_max[_qp];
                _kappa[_qp] = _kappa_max[_qp];
                _beta[_qp] = _beta_max[_qp];
            }
            else {
                _beta[_qp] = std::min(_beta_max[_qp], _Ce*_h*_h*(std::fabs(_Dalpha)+_jump_alpha) / _norm_alpha);
                _kappa[_qp] = std::min( _kappa_max[_qp], _kappa_e );
                _mu[_qp] = std::min( _kappa_max[_qp], _mu_e );
            }
            break;
        default:
            mooseError("The viscosity type entered in the input file is not implemented.");
            break;
    }
}
