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
    // Bool for viscosity coefficient of void fraction equation:
    params.addParam<bool>("isShock", false, "employ definition of the viscosity coefficient for weak shock.");
    params.addParam<bool>("areViscEqual", false, "the viscosity coefficients are equal in the shock region.");
    // Aux variables:
    params.addRequiredCoupledVar("velocity_x", "x component of the velocity");
    params.addCoupledVar("velocity_y", "y component of the velocity");
    params.addCoupledVar("velocity_z", "z component of the velocity");
    params.addRequiredCoupledVar("pressure", "pressure of the fluid");
    params.addRequiredCoupledVar("density", "density of the fluid: rho");
    params.addCoupledVar("vf_liquid", "liquid void fraction.");
    params.addCoupledVar("area", 1., "area, cross-section.");
    // Jumps:
    params.addCoupledVar("jump_grad_press", "jump of pressure gradient");
    params.addCoupledVar("jump_grad_dens", "jump of density gradient");
    params.addCoupledVar("jump_grad_alpha", "jump of alpha gradient");
    // Constant parameter:
    params.addParam<double>("Cmax", 0.5, "Coefficient for first-order viscosity");
    params.addParam<double>("Ce", 1., "Coefficient for residual");
    params.addParam<double>("Cjump_liquid", 1., "Coefficient for jump for liquid phase");
    params.addParam<double>("Cjump_gas", 1., "Coefficient for jump for gas phase");
    params.addParam<double>("Calpha", 1., "Coefficient for alpha");
    // Coefficients for 'sigma' function:
    params.addParam<Real>("a_coeff", 4., "Coefficient for sigma function");
    params.addParam<Real>("Mthres", 0.005, "Threshold Mach number");
    // Userobject:
    params.addRequiredParam<UserObjectName>("eos", "Equation of state");
    // PPS names:
    params.addParam<std::string>("rhov2_PPS_name", "name of the pps computing rho*vel*vel");
    params.addRequiredParam<std::string>("alpha_PPS_name", "name of the pps for alpha");
    return params;
}

ComputeViscCoeff::ComputeViscCoeff(const std::string & name, InputParameters parameters) :
    Material(name, parameters),
    // Declare viscosity types
    _visc_name(getParam<std::string>("viscosity_name")),
    _visc_type("LAPIDUS, FIRST_ORDER, FIRST_ORDER_MACH, ENTROPY, INVALID", _visc_name),
    // Function Mach number:
    _fct_of_mach_name(getParam<std::string>("function_of_mach")),
    _fct_of_mach_type("MACH, SQRT_MACH, FCT_OF_MACH, INVALID", _fct_of_mach_name),
    // Boolean for phase:
    _isLiquid(getParam<bool>("isLiquid")),
    _isJumpOn(getParam<bool>("isJumpOn")),
    // Bool for viscosity coefficient of void fraction equation:
    _isShock(getParam<bool>("isShock")),
    _areViscEqual(getParam<bool>("areViscEqual")),
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
    // Get interfacial velocity
    _velI(getMaterialProperty<RealVectorValue>("interfacial_velocity")),
    // Get parameters: Cmax, Ce, Cjump:
    _Cmax(getParam<double>("Cmax")),
    _Ce(getParam<double>("Ce")),
    _Cjump(_isLiquid ? getParam<double>("Cjump_liquid") : getParam<double>("Cjump_gas")),
    _Calpha(getParam<double>("Calpha")),
    // Coefficients for 'sigma' function:
    _a_coeff(getParam<Real>("a_coeff")),
    _Mthres(getParam<Real>("Mthres")),
    // UserObject:
    _eos(getUserObject<EquationOfState>("eos")),
    // PPS name:
    _rhov2_pps_name(getParam<std::string>("rhov2_PPS_name")),
    _alpha_pps_name(getParam<std::string>("alpha_PPS_name"))
{
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
    
    /** Compute the first order viscosity and the mach number: **/
    Real c = std::sqrt(_eos.c2_from_p_rho(_rho[_qp], _pressure[_qp]));
    RealVectorValue vel(_vel_x[_qp], _vel_y[_qp], _vel_z[_qp]);
    _mu_max[_qp] = _Cmax*h*vel.size();
    _kappa_max[_qp] = _Cmax*h*(vel.size() + c);
    _beta_max[_qp] = _Cmax*h*_velI[_qp].size();
    
    /** Compute the Mach number: **/
    Real Mach = std::min(vel.size()/c, 1.);
//    Real Mach2 = _areViscEqual ? Mach*Mach: 1.;
    Real Mach2 = Mach*Mach;
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
    
    // Postprocessors:
    Real rhov2_pps = std::max(getPostprocessorValueByName(_rhov2_pps_name), eps);
    Real alpha_var = getPostprocessorValueByName(_alpha_pps_name);
    
    // Initialyze some variables used in the switch statement:
//    h = _current_elem->hmin()/_qrule->get_order();
    Real weight0, weight1, weight2;
    Real mu_e, kappa_e, beta_e;
    Real jump, residual, norm, sigma;
    
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
            _mu[_qp] = _kappa_max[_qp];
            _kappa[_qp] = _kappa_max[_qp];
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
            if (_isJumpOn)
                jump = _Calpha*std::fabs(_velI[_qp].size()*_jump_grad_alpha[_qp]);
            else
                jump = _Calpha*std::fabs(_velI[_qp].size()*_grad_alpha_l[_qp](0));
//            norm = std::min(_alpha_l[_qp], std::fabs(1.-_alpha_l[_qp]));
            norm = alpha_var;
            beta_e = h*h*(std::fabs(residual)+jump) / norm;
//            beta_e += h*h*std::fabs(vel*_grad_area[_qp])/_area[_qp];
            
        /** Compute viscosity coefficient for continuity, momentum and energy equations: **/
            // Entropy residual:
            residual = 0.;
            residual = vel*_grad_press[_qp];
            residual += (weight0*_pressure[_qp]+weight1*_pressure_old[_qp]+weight2*_pressure_older[_qp]);
            residual -= c*c*vel*_grad_rho[_qp];
            residual -= c*c*(weight0*_rho[_qp]+weight1*_rho_old[_qp]+weight2*_rho_older[_qp]);
            residual *= _Ce;
            
            if (_isShock) // non-isentropic flow.
            {
                // Compute the jumps for mu_e:
                if (_isJumpOn)
                    jump = _Cjump*vel.size()*std::max( _jump_grad_press[_qp], Mach2*c*c*_jump_grad_dens[_qp] );
                else
                    jump = _Cjump*vel.size()*std::max( _grad_press[_qp].size(), Mach2*c*c*_grad_rho[_qp].size() );
                
                // Compute high-order viscosity coefficient mu_e:
                norm = 0.5 * std::max( (1.-Mach)*rhov2_pps, _rho[_qp]*std::min(vel.size_sq(), c*c) );
                mu_e = h*h*(std::fabs(residual) + jump) / norm;
                mu_e += h*h*_pressure[_qp] * std::fabs(vel * _grad_area[_qp]) / ( _area[_qp] * norm );
                                
                // Compute the jumps for kappa_e:
                if (_isJumpOn)
                    jump = _Cjump*vel.size()*std::max( _jump_grad_press[_qp], c*c*_jump_grad_dens[_qp] );
                else
                    jump = _Cjump*vel.size()*std::max( _grad_press[_qp].size(), c*c*_grad_rho[_qp].size() );

                // Compute high-order viscosity coefficient kappa_e:
//                norm = 0.5*( std::fabs(1.-Mach)*_rho[_qp]*c*c + Mach*_rho[_qp]*std::min(vel.size_sq(), c*c) );
                if (!_areViscEqual)
                    norm = 0.5*_rho[_qp]*c*c;
                kappa_e = h*h*(std::fabs(residual) + jump) / norm;
                kappa_e += h*h*_pressure[_qp] * std::fabs(vel * _grad_area[_qp]) / ( _area[_qp] * norm );
            }
            else // with function sigma.
            {
                // Compute the jumps:
                if (_isJumpOn)
                    jump = _Cjump*vel.size()*std::max( _jump_grad_press[_qp], c*c*_jump_grad_dens[_qp] );
                else
                    jump = _Cjump*vel.size()*std::max( _grad_press[_qp].size(), c*c*_grad_rho[_qp].size() );
                
                // Compute the function sigma:
                sigma = 0.5 * std::tanh(_a_coeff*(Mach-_Mthres));
                sigma += 0.5 * std::abs(std::tanh(_a_coeff*(Mach-_Mthres)));
                
                // Compute mu_e:
                norm = (1.-sigma) * _rho[_qp] * c * c;
                norm += sigma * _rho[_qp] * vel.size_sq();
                norm *= 0.5;
                mu_e = h*h*(std::fabs(residual) + jump) / norm;
                
                // Compute kappa_e:
                norm = 0.5 * _rho[_qp] * c * c;
                kappa_e = h*h*(std::fabs(residual) + jump) / norm;
                
//                // Compute high-order viscosity coefficients:
//                norm = 0.5*( std::fabs(1.-Mach)*_rho[_qp]*c*c + Mach*_rho[_qp]*vel.size_sq() );
//                kappa_e = h*h*std::max(std::fabs(residual), jump) / norm;
////                kappa_e += h*h*_pressure[_qp] * std::fabs(vel * _grad_area[_qp]) / ( _area[_qp] * norm );
//                kappa_e += h*h*std::fabs(vel*_grad_area[_qp])/_area[_qp];
//                mu_e = kappa_e;
            }
            
        /** Compute the viscosity coefficients: **/
            if (_t_step == 1)
            {
                _mu[_qp] = 2*_kappa_max[_qp];
                _kappa[_qp] = 2*_kappa_max[_qp];
                _beta[_qp] = 2*_beta_max[_qp];
            }
            else
            {
                _beta[_qp] = std::min(_beta_max[_qp], beta_e);
                _kappa[_qp] = std::min( _kappa_max[_qp], kappa_e );
                _mu[_qp] = std::min( _kappa_max[_qp], mu_e );
            }
            break;
        default:
            mooseError("The viscosity type entered in the input file is not implemented.");
            break;
    }
}
