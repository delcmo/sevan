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
    params.addParam<bool>("isJumpOn", true, "use jump or gradients.");
    params.addParam<bool>("isVariableArea", false, "Smooth variable area");
    // Bool for viscosity coefficient of void fraction equation:
    params.addParam<bool>("useLiqViscForVF", false, "set beta equal to mu for liquid phase");
    // Aux variables:
    params.addRequiredCoupledVar("velocity_x", "x component of the velocity");
    params.addCoupledVar("velocity_y", "y component of the velocity");
    params.addCoupledVar("velocity_z", "z component of the velocity");
    params.addRequiredCoupledVar("pressure", "pressure of the fluid");
    params.addRequiredCoupledVar("density", "density of the fluid: rho");
    params.addRequiredCoupledVar("internal_energy", "internal energy");
    params.addCoupledVar("jump_grad_press", "jump of pressure gradient");
    params.addCoupledVar("jump_grad_dens", "jump of density gradient");
    params.addCoupledVar("jump_grad_alpha", "jump of alpha gradient");
    params.addRequiredCoupledVar("norm_velocity", "norm of the velocity vector");
    params.addCoupledVar("vf_liquid", "liquid void fraction.");
    params.addRequiredCoupledVar("area", "area, cross-section.");
    // Constant parameter:
    params.addParam<double>("Ce", 1., "Coefficient for residual");
    params.addParam<double>("Cjump", 1., "Coefficient for jump");
    params.addParam<double>("Calpha", 1., "Coefficient for alpha");
    // Userobject:
    params.addRequiredParam<UserObjectName>("eos", "Equation of state");
    // PPS names:
    params.addParam<std::string>("rhov2_PPS_name", "name of the pps computing rho*vel*vel");
    params.addParam<std::string>("rhoc2_PPS_name", "name of the pps computing rho*c*c");
    params.addRequiredParam<std::string>("alpha_PPS_name", "name of the pps for alpha");
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
    _isJumpOn(getParam<bool>("isJumpOn")),
    _isVariableArea(getParam<bool>("isVariableArea")),
    // Bool for viscosity coefficient of void fraction equation:
    _useLiqViscForVF(getParam<bool>("useLiqViscForVF")),
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
    // Internal energy:
    _rhoe(coupledValue("internal_energy")),
    // Jump of pressure, density and alpha gradients:
    _jump_grad_press(isCoupled("jump_grad_press") ? coupledValue("jump_grad_press") : _zero),
    _jump_grad_dens(isCoupled("jump_grad_dens") ? coupledValue("jump_grad_dens") : _zero),
    _jump_grad_alpha(isCoupled("jump_grad_alpha") ? coupledValue("jump_grad_alpha") : _zero),
    // Norm of velocity vector:
    _norm_vel(coupledValue("norm_velocity")),
    // Area
    _area(coupledValue("area")),
    _grad_area(coupledGradient("area")),
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
    _rhov2_pps_name(getParam<std::string>("rhov2_PPS_name")),
    _rhoc2_pps_name(getParam<std::string>("rhoc2_PPS_name")),
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
    Real h = _current_elem->hmin();
    Real eps = std::sqrt(std::numeric_limits<Real>::min());
    
    // Compute the first order viscosity and the mach number:
    Real c = std::sqrt(_eos.c2_from_p_rho(_rho[_qp], _pressure[_qp]));
    _mu_max[_qp] = 0.5*h*_norm_vel[_qp];
    _kappa_max[_qp] = 0.5*h*(_norm_vel[_qp] + c);
    _beta_max[_qp] = 0.5*h*_velI[_qp].size();
    Real Mach = std::min(_norm_vel[_qp]/c, 1.);
    Real fct_of_mach = Mach;
    switch (_fct_of_mach_type) {
        case MACH:
            fct_of_mach = std::min(Mach, 1.);
            break;
        case SQRT_MACH:
            fct_of_mach = std::min(std::sqrt(Mach), 1.);
            break;
        case FCT_OF_MACH:
            fct_of_mach = std::min(Mach*std::sqrt(4+(1.-Mach*Mach)*(1.-Mach*Mach)) / (1.+Mach*Mach),1.);
            break;
        default:
            mooseError("The function with name: \"" << _fct_of_mach_name << "\" is not supported in the \"ComputeViscCoeff\" type of material.");
    }
    
    // Postprocessors
    Real rhov2_pps = std::max(getPostprocessorValueByName(_rhov2_pps_name), eps);
    Real rhocv_pps = std::max(getPostprocessorValueByName(_rhocv_pps_name), eps);
    Real rhoc2_pps = std::max(getPostprocessorValueByName(_rhoc2_pps_name), eps);
    Real press_pps = std::max(getPostprocessorValueByName(_press_pps_name), eps);
    Real alpha_var = getPostprocessorValueByName(_alpha_pps_name);
    
    // Initialyze some variables used in the switch statement:
    Real weight0, weight1, weight2;
    Real jump, kappa_e1, kappa_e2;
    Real kappa_e, beta_e, residual, norm;
    RealVectorValue vel(_vel_x[_qp], _vel_y[_qp], _vel_z[_qp]);
    
    // Switch statement over viscosity type:
    switch (_visc_type) {
        case LAPIDUS:
            if (_t_step == 1) {
                _mu[_qp] = _mu_max[_qp];
                _kappa[_qp] = _kappa_max[_qp];
                _beta[_qp] = _beta_max[_qp];
            }
            else {
                _mu[_qp] = _Ce*h*h*std::fabs(_grad_vel_x[_qp](0));
                _kappa[_qp] = _mu[_qp];
                _beta[_qp] = _mu[_qp];
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
            _mu[_qp] = Mach*Mach*_mu_max[_qp];
            _kappa[_qp] = _kappa_max[_qp];
            _beta[_qp] = _beta_max[_qp];
            break;
        case ENTROPY:
            // Compute the weights for BDF2
            if (_t_step > 1)
            {
                weight0 = (2.*_dt+_dt_old)/(_dt*(_dt+_dt_old));
                weight1 = -(_dt+_dt_old)/(_dt*_dt_old);
                weight2 = _dt/(_dt_old*(_dt+_dt_old));
            }
            else
            {
                weight0 =  1. / _dt;
                weight1 = -1. / _dt;
                weight2 = 0.;
            }
            
        /** Compute viscosity coefficient for void fraction equation: **/
            residual = 0.;
            residual = _velI[_qp]*_grad_alpha_l[_qp];
            residual += (weight0*_alpha_l[_qp]+weight1*_alpha_l_old[_qp]+weight2*_alpha_l_older[_qp]);
            residual *= _Ce;
            jump = _Calpha*std::fabs(_velI[_qp].size()*_grad_alpha_l[_qp](0));
            norm = alpha_var;
            beta_e = h*h*(std::fabs(residual)+jump) / norm;
            
        /** Compute viscosity coefficient for continuity, momentum and energy equations: **/
            // Entropy residual:
            residual = 0.;
            residual = vel*_grad_press[_qp];
            residual += (weight0*_pressure[_qp]+weight1*_pressure_old[_qp]+weight2*_pressure_older[_qp]);
            residual -= c*c*vel*_grad_rho[_qp];
            residual -= c*c*(weight0*_rho[_qp]+weight1*_rho_old[_qp]+weight2*_rho_older[_qp]);
            residual *= _Ce;
            
            // Compute the jumps:
            if (_isJumpOn)
                jump = _Cjump*_norm_vel[_qp]*std::max( _jump_grad_press[_qp], c*c*_jump_grad_dens[_qp] );
            else
                jump = _Cjump*_norm_vel[_qp]*std::max( _grad_press[_qp].size(), c*c*_grad_rho[_qp].size() );
            
            // Compute kappa1 and kappa2:
//            norm = std::max((1.-Mach)*rhov2_pps, 0.5*_rho[_qp]*std::min(_norm_vel[_qp]*_norm_vel[_qp], c*c));
//            std::cout<<"residual="<<_dt*std::fabs(residual)/_pressure[_qp]<<" and "<<Mach<<std::endl;
//            if (_isVariableArea)
//                norm = std::max(_rho[_qp]*std::min(_norm_vel[_qp]*_norm_vel[_qp], c*c), (1.-Mach)*rhov2_pps );
//            else if (_dt*std::fabs(residual)/_pressure[_qp] > Mach)
//                norm = 0.5 * _rho[_qp] * std::min(_norm_vel[_qp]*_norm_vel[_qp], c*c);
//            else {
//                norm = std::fabs(1.-Mach) * _rho[_qp] * c*c;
//                norm += Mach * _rho[_qp] * std::min(_norm_vel[_qp]*_norm_vel[_qp], c*c);
//                norm *= 0.5;
//            }
//            norm = std::max(_rho[_qp]*std::min(_norm_vel[_qp]*_norm_vel[_qp], c*c), (1.-Mach)*rhov2_pps );
            norm = 0.5*( std::fabs(1.-Mach)*_rho[_qp]*c*c + Mach*_rho[_qp]*_norm_vel[_qp]*_norm_vel[_qp] );
//            norm = 0.5*_rho[_qp]*c*c;
            kappa_e1 = h*h*std::max(std::fabs(residual), jump) / norm;
//            norm = 0.5*_rho[_qp]*c*c;
            kappa_e2 = h*h*std::max(std::fabs(residual), jump ) / norm;
            kappa_e2 += h*h*_pressure[_qp]*_norm_vel[_qp]*std::fabs(_grad_area[_qp].size())/_area[_qp] / norm;
            
        /** Compute the viscosity coefficients: **/
            if (_t_step == 1)
            {
                _mu[_qp] = _mu_max[_qp];
                _kappa[_qp] = _kappa_max[_qp];
                _beta[_qp] = _beta_max[_qp];
            }
            else
            {
                _beta[_qp] = std::min(_beta_max[_qp], beta_e);
                _kappa[_qp] = std::min( _kappa_max[_qp], kappa_e1 );
                _mu[_qp] = std::min( _kappa_max[_qp], kappa_e2 );
            }
            break;
        default:
            mooseError("The viscosity type entered in the input file is not implemented.");
            break;
    }
}
