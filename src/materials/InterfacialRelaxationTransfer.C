#include "InterfacialRelaxationTransfer.h"
/* This function computes the interfacial Relaxation (PI, velI, PI_bar and velI_bar) and the relaxation parameters (mu and lambda) for the 7 equations model. 
 It also computes the wall heat transfer and friction parameters for each phase.*/
template<>
InputParameters validParams<InterfacialRelaxationTransfer>()
{
  InputParameters params = validParams<Material>();
    params.addParam<std::string>("inter_def_name", "BERRY", "Choose definition to compute interfacial variables.");
    params.addParam<Real>("xi_ambrosso", 0.5, "value for definition of interfacial variables");
    // Boolean for mass and heat transfers:
    params.addParam<bool>("isMassOn", false, "are the mass transfer terms are on?");
    params.addParam<bool>("isHeatOn", false, "are the heat transfer terms are on?");
    params.addParam<bool>("isWallHeatOn", false, "are the wall heat transfer terms are on?");
    params.addParam<bool>("isWallFrictOn", false, "are the wall friction terms are on?");
    // Boolean for interfacial relaxation parameters:
    params.addParam<bool>("isPressRelOn", true, "is pressure relaxation on?");
    params.addParam<bool>("isVelRelOn", true, "is velocity relaxation on?");
    // Aux Relaxation for liquid phase:
    params.addRequiredCoupledVar("velocity_x_liq", "x component of the liquid velocity");
    params.addCoupledVar("velocity_y_liq", "y component of the liquid velocity");
    params.addCoupledVar("velocity_z_liq", "z component of the liquid velocity");
    params.addRequiredCoupledVar("pressure_liq", "pressure of the liquid");
    params.addRequiredCoupledVar("density_liq", "density of the liquid: rho");
    params.addRequiredCoupledVar("vf_liquid", "LIQUID void fraction");
    // Aux Relaxation for gas phase:
    params.addRequiredCoupledVar("velocity_x_gas", "x component of the gas velocity");
    params.addCoupledVar("velocity_y_gas", "y component of the gas velocity");
    params.addCoupledVar("velocity_z_gas", "z component of the gas velocity");
    params.addRequiredCoupledVar("pressure_gas", "pressure of the gas");
    params.addRequiredCoupledVar("density_gas", "density of the gas: rho");
    // Constant:
    params.addParam<Real>("Aint", 0., "Specific interfacial area");
    params.addParam<Real>("wall_heat_liq_value", 0., "Wall heat transfer value for liquid phase.");
    params.addParam<Real>("wall_heat_gas_value", 0., "Wall heat transfer value for gas phase.");
    params.addParam<Real>("wall_frict_liq_value", 0., "Wall friction value for liquid phase.");
    params.addParam<Real>("wall_frict_gas_value", 0., "Wall friction value for gas phase.");
    params.addParam<Real>("Twall", 0., "Constant wall temperature.");
    // Userobject:
    params.addRequiredParam<UserObjectName>("eos_liq", "Liquid equation of state");
    params.addRequiredParam<UserObjectName>("eos_gas", "Gas equation of state");
    return params;
}

InterfacialRelaxationTransfer::InterfacialRelaxationTransfer(const std::string & name, InputParameters parameters) :
    Material(name, parameters),
    // Definition for interfacial variables:
    _interf_def_name(getParam<std::string>("inter_def_name")),
    _interf_def("BERRY, AMBROSSO, LIANG", _interf_def_name),
    _xi(getParam<Real>("xi_ambrosso")),
    // Boolean for mass and heat transfers:
    _isMassOn(getParam<bool>("isMassOn")),
    _isHeatOn(getParam<bool>("isHeatOn")),
    _isWallHeatOn(getParam<bool>("isWallHeatOn")),
    _isWallFrictOn(getParam<bool>("isWallFrictOn")),
    // Boolean for interfacial relaxation parameters:
    _isPressRelOn(getParam<bool>("isPressRelOn")),
    _isVelRelOn(getParam<bool>("isVelRelOn")),
    // Aux variable for liquid phase:
    _vel_x_l(coupledValue("velocity_x_liq")),
    _vel_y_l(_mesh.dimension()>=2 ? coupledValue("velocity_y_liq") : _zero),
    _vel_z_l(_mesh.dimension()==3 ? coupledValue("velocity_z_liq") : _zero),
    _pressure_l(coupledValue("pressure_liq")),
    _rho_l(coupledValue("density_liq")),
    _alpha_l(coupledValue("vf_liquid")),
    _grad_alpha_l(coupledGradient("vf_liquid")),
    // Aux variable for gas phase:
    _vel_x_g(coupledValue("velocity_x_gas")),
    _vel_y_g(_mesh.dimension()>=2 ? coupledValue("velocity_y_gas") : _zero),
    _vel_z_g(_mesh.dimension()==3 ? coupledValue("velocity_z_gas") : _zero),
    _pressure_g(coupledValue("pressure_gas")),
    _rho_g(coupledValue("density_gas")),
    // Declare interfacial variables:
    _Aint(declareProperty<Real>("interfacial_area")),
    _PI(declareProperty<Real>("interfacial_pressure")),
    _velI(declareProperty<RealVectorValue>("interfacial_velocity")),
    _velI_norm(declareProperty<Real>("interfacial_velocity_norm")),
    _PI_bar(declareProperty<Real>("average_interfacial_pressure")),
    _velI_bar(declareProperty<RealVectorValue>("average_interfacial_velocity")),
    _tempI(declareProperty<Real>("interfacial_temperature")),
    _rhoI(declareProperty<Real>("interfacial_density")),
    _EI_liq(declareProperty<Real>("liquid_interfacial_energy")),
    _EI_gas(declareProperty<Real>("gas_interfacial_energy")),
    // Declare relaxation parameters:
    _P_rel(declareProperty<Real>("pressure_relaxation")),
    _vel_rel(declareProperty<Real>("velocity_relaxation")),
    // Heat transfer coefficient for liquid and gas phases:
    _ht_liq(declareProperty<Real>("liquid_heat_transfer")),
    _ht_gas(declareProperty<Real>("gas_heat_transfer")),
    // Mass transfer coefficient:
    _Omega_gas(declareProperty<Real>("mass_transfer")),
    // Constant:
    _Aint_max(getParam<Real>("Aint")),
    _wall_heat_liq_value(getParam<Real>("wall_heat_liq_value")),
    _wall_heat_gas_value(getParam<Real>("wall_heat_gas_value")),
    _wall_frict_liq_value(getParam<Real>("wall_frict_liq_value")),
    _wall_frict_gas_value(getParam<Real>("wall_frict_gas_value")),
    _Twall(getParam<Real>("Twall")),
    // UserObject:
    _eos_liq(getUserObject<EquationOfState>("eos_liq")),
    _eos_gas(getUserObject<EquationOfState>("eos_gas")),
    // Wall heat transfer for liquid and gas phases:
    _wall_ht_liq(declareProperty<Real>("wall_heat_transfer_liq")),
    _wall_ht_gas(declareProperty<Real>("wall_heat_transfer_gas")),
    // Wall friction parameters for liquid and gas phases:
    _wall_frict_liq(declareProperty<Real>("wall_friction_liq")),
    _wall_frict_gas(declareProperty<Real>("wall_friction_gas")),
    // Interfacial friction parameters for liquid and gas phases:
    _interf_frict_liq(declareProperty<Real>("interfacial_friction_liq")),
    _interf_frict_gas(declareProperty<Real>("interfacial_friction_gas")),
    // Wall temperature:
    _wall_temp(declareProperty<Real>("wall_temperature"))
{
}

void
InterfacialRelaxationTransfer::computeQpProperties()
{
    // Initialize velocity vectors for each phase:
    RealVectorValue _vel_l(_vel_x_l[_qp], _vel_y_l[_qp], _vel_z_l[_qp]);
    RealVectorValue _vel_g(_vel_x_g[_qp], _vel_y_g[_qp], _vel_z_g[_qp]);

    /*******************************************************************************/
    /********** Compute interfacial variables and relaxation coefficients: *********/
    /*******************************************************************************/
    // Compute the speed of sound for each phase:
    Real _c2_l = _eos_liq.c2_from_p_rho(_rho_l[_qp], _pressure_l[_qp]);
    Real _c2_g = _eos_gas.c2_from_p_rho(_rho_g[_qp], _pressure_g[_qp]);
    
    // Compute the impedences for each phase:
    Real _Z_l = _rho_l[_qp] * std::sqrt(_c2_l);
    Real _Z_g = _rho_g[_qp] * std::sqrt(_c2_g);
    Real _sum_Z = _Z_l + _Z_g;
    
    // Compute the interfacial variables according to the definition chose in input file:
    Real beta, mu, temp_l, temp_g;
    RealVectorValue _n(_grad_alpha_l[_qp](0), _grad_alpha_l[_qp](1), _grad_alpha_l[_qp](2));
    switch (_interf_def)
    {
        case BERRY:
            // Compute unit vector based on gradient of liquid void fraction:
            if ( _mesh.dimension() == 1 )
            {
                _n(0) = _grad_alpha_l[_qp](0) >= 0 ? 1. : -1.;
                _n(1) = 0.; _n(2) = 0.;
            }
            else
            {
                Real _eps = std::sqrt(std::numeric_limits<Real>::min());
                if (_n.size() <= 1e-4) {
                    _n(0) = 0.; _n(1) = 0.; _n(2) = 0.;
                }
                else {
                    _n = _n / (_n.size() + _eps); }
            }
            
            if (_Aint_max == 0)
            {
                // Compute the average interfacial and relaxation parameters:
                _PI_bar[_qp] = ( _Z_g*_pressure_l[_qp] + _Z_l*_pressure_g[_qp] ) / _sum_Z;
                _velI_bar[_qp] = ( _Z_l*_vel_l + _Z_g*_vel_g ) / _sum_Z;
                
                // Compute interfacial Relaxation parameters:
                _PI[_qp] = _PI_bar[_qp] + _Z_l*_Z_g/_sum_Z * _n*(_vel_g-_vel_l);
                _velI[_qp] = _velI_bar[_qp] + _n * (_pressure_g[_qp]-_pressure_l[_qp])/_sum_Z;
                
                _velI[_qp] = _velI_bar[_qp];
                _Aint[_qp] = 0.; _P_rel[_qp] = 0.; _vel_rel[_qp] = 0.;
            }
            else
            {
                // Compute the average interfacial Relaxation parameters:
                _PI_bar[_qp] = ( _Z_g*_pressure_l[_qp] + _Z_l*_pressure_g[_qp] ) / _sum_Z;
                _velI_bar[_qp] = ( _Z_l*_vel_l + _Z_g*_vel_g ) / _sum_Z;
                
                // Compute interfacial Relaxation parameters:
                _PI[_qp] = _PI_bar[_qp] + _Z_l*_Z_g/_sum_Z * _n*(_vel_g-_vel_l);
                _velI[_qp] = _velI_bar[_qp] + _n * (_pressure_g[_qp]-_pressure_l[_qp])/_sum_Z;
            }
            break;
        case AMBROSSO:
            // Compute interfacial velocity:
            beta = _xi*_alpha_l[_qp]*_rho_l[_qp];
            beta *= 1./(_xi*_alpha_l[_qp]*_rho_l[_qp] + (1.-_xi)*(1.-_alpha_l[_qp])*_rho_g[_qp]);
            _velI_bar[_qp] = beta*_vel_l + (1.-beta)*_vel_g;
            _velI[_qp] = _velI_bar[_qp];
            // Compute interfacial pressure:
            temp_l = _eos_liq.temperature_from_p_rho(_pressure_l[_qp], _rho_l[_qp]);
            temp_g = _eos_gas.temperature_from_p_rho(_pressure_g[_qp], _rho_g[_qp]);
            mu = (1.-beta)*temp_g/(beta*temp_l+(1.-beta)*temp_g);
            _PI_bar[_qp] = mu*_pressure_l[_qp] + (1.-mu)*_pressure_g[_qp];
            _PI[_qp] = _PI_bar[_qp];
            break;
        case LIANG:
            // Compute intefacial velocity:
            _velI_bar[_qp] = _alpha_l[_qp]*_rho_l[_qp]*_vel_l + (1.-_alpha_l[_qp])*_rho_g[_qp]*_vel_g;
            _velI_bar[_qp] *= 1./(_alpha_l[_qp]*_rho_l[_qp] + (1.-_alpha_l[_qp])*_rho_g[_qp]);
            _velI[_qp] = _velI_bar[_qp];
            // Compute interfacial pressure:
            _PI_bar[_qp] = _alpha_l[_qp]*_pressure_l[_qp] + (1.-_alpha_l[_qp])*_pressure_g[_qp];
            break;
            
        default:
            break;
    }

    // Compute the relaxation parameters:
    _Aint[_qp] = _Aint_max;//*(6.75*(1-_alpha_l[_qp])*(1-_alpha_l[_qp])*_alpha_l[_qp]);
    _P_rel[_qp] = _isPressRelOn ? _Aint[_qp]/_sum_Z : 0.; /*(mu)*/
    _vel_rel[_qp] = _isVelRelOn ? 0.5*_Aint[_qp]*_Z_g*_Z_l/_sum_Z : 0.; /*(lambda)*/
    
    // Compute the norm of the interfacial velocity for output only
    _velI_norm[_qp] = _velI[_qp].size();

    /***************************************************/
    /***** Newton solve for computing TI from PI: ******/
    /***************************************************/
    
    Real _temp = _eos_liq.temperature_from_p_rho(_pressure_l[_qp], _rho_l[_qp]);
    if (_isMassOn == true) {
        // Define some parameters used in the local Newton solve: (DEM paper)
        Real _A = (_eos_liq.Cp() - _eos_gas.Cp() + _eos_gas.qcoeff_prime() - _eos_liq.qcoeff_prime()) / (_eos_gas.Cp() - _eos_gas.Cv());
        Real _B = (_eos_liq.qcoeff() - _eos_gas.qcoeff()) / (_eos_gas.Cp() - _eos_gas.Cv());
        Real _C = (_eos_gas.Cp() - _eos_liq.Cp()) / (_eos_gas.Cp() - _eos_gas.Cv());
        Real _D = (_eos_liq.Cp() - _eos_liq.Cv()) / (_eos_gas.Cp() - _eos_gas.Cv());

        // Compute the constant residual _R:
        Real _p_term_liq = std::log(_PI[_qp]+_eos_liq.Pinf());
        Real _p_term_gas = std::log(_PI[_qp]+_eos_gas.Pinf());
        Real _R = _A + _D * _p_term_liq - _p_term_gas;

        // Initialyze some values and pick a guess for the temperature:
        Real _f_norm = 1;
        Real _f = 0.0;
        Real _f_prime = 0.0;

        // Newton solve:
        while ( std::fabs(_f_norm) > 1e-5)
        {
            _f = _R + _B / _temp + _C * std::log(_temp);
            _f_prime = _C / _temp - _B / (_temp*_temp);
            _temp = _temp - _f / _f_prime;
            _f_norm = _f / _f_prime;
        }
    }
    _tempI[_qp] = _temp;

    // Compute the interfacial density:
    //    _rhoI[_qp] = _eos_liq.rho_from_p_T(_pressure_l[_qp], _tempI[_qp]);
    _rhoI[_qp] = _eos_liq.rho_from_p_T(_PI[_qp], _tempI[_qp]);

    /*********************************************************/
    /** Compute the interfacial heat transfer coefficients: **/
    /*********************************************************/
    
    Real _Lv = (_eos_gas.Cp()-_eos_liq.Cp())*_tempI[_qp] + (_eos_gas.qcoeff()-_eos_liq.qcoeff());

    // Compute the interfacial energy for liquid and gas phases:
    _EI_liq[_qp] = _eos_liq.Cp()*_tempI[_qp] + _eos_liq.qcoeff() + 0.5*_velI[_qp]*_velI[_qp];
    _EI_gas[_qp] = _eos_gas.Cp()*_tempI[_qp] + _eos_gas.qcoeff() + 0.5*_velI[_qp]*_velI[_qp];

    // Compute the heat transfer coefficient for liquid phase:
    Real _radius = 1.;
    if (_Aint_max != 0)
        _radius = 3*(1-_alpha_l[_qp])/_Aint[_qp];

    _ht_liq[_qp] = _isHeatOn ? 5*_eos_liq.k()/_radius : 0.;

    // Compute the heat transfer coefficient for gas phase:
    Real _L = 2*_radius;
    Real _Re = _rho_g[_qp]*_L*std::fabs(_vel_l.size()-_vel_g.size())/_eos_gas.visc();
    Real _Pr = _eos_gas.visc()*_eos_gas.Cp()/_eos_gas.k();
    Real _Nu = 2 + 0.6*std::pow(_Re,0.5) * std::pow(_Pr, 0.33);

    _ht_gas[_qp] = _isHeatOn ? _eos_gas.k()*_Nu/_L : 0.;

    // Compute the mass transfer coefficient between phases:
    Real _temp_liq = _eos_liq.temperature_from_p_rho(_pressure_l[_qp], _rho_l[_qp]);
    Real _temp_gas = _eos_gas.temperature_from_p_rho(_pressure_g[_qp], _rho_g[_qp]);

    _Omega_gas[_qp] = _isMassOn ? ( _ht_liq[_qp]*(_temp_liq-_tempI[_qp])+_ht_gas[_qp]*(_temp_gas-_tempI[_qp]) ) / _Lv : 0.;

    /****************************************************************/
    /* Compute the wall heat transfer and friction for liquid phase */
    /****************************************************************/
    
    _wall_ht_liq[_qp] = 0.;
    _wall_frict_liq[_qp] = 0.;
    if (_isWallHeatOn)
        _wall_ht_liq[_qp] = _wall_heat_liq_value;
    if (_isWallFrictOn)    
        _wall_frict_liq[_qp] = _wall_frict_liq_value;

    /*************************************************************/
    /* Compute the wall heat transfer and friction for gas phase */
    /*************************************************************/
    
    _wall_ht_gas[_qp] = 0.;
    _wall_frict_gas[_qp] = 0.;
    if (_isWallHeatOn)
        _wall_ht_gas[_qp] = _wall_heat_gas_value;
    if (_isWallFrictOn)
        _wall_frict_gas[_qp] = _wall_frict_gas_value;
    
    /*****************************************************************/
    /* Compute the interfacial friction coefficient for liquid phase */
    /*****************************************************************/
    
    _interf_frict_liq[_qp] = 0.;
    
    /**************************************************************/
    /* Compute the interfacial friction coefficient for gas phase */
    /**************************************************************/

    _interf_frict_gas[_qp] = 0.;
    
    /********************************/
    /* Compute the wall temperature */
    /********************************/
    
    _wall_temp[_qp] = 0.;
    if (_isWallHeatOn) {
        _wall_temp[_qp] = _Twall;
    }
}