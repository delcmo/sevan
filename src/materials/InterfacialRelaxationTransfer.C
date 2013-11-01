#include "InterfacialRelaxationTransfer.h"
/* This function computes the interfacial Relaxation (PI, velI, PI_bar and velI_bar) and the relaxation parameters (mu and lambda) for the 7 equations model. */
template<>
InputParameters validParams<InterfacialRelaxationTransfer>()
{
  InputParameters params = validParams<Material>();
    // Boolean for mass and heat transfers:
    params.addParam<bool>("isMassOn", false, "are the mass transfer terms are on?");
    params.addParam<bool>("isHeatOn", false, "are the heat transfer terms are on?");
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
    params.addParam<Real>("Aint", 0, "Specific interfacial area");
    // Userobject:
    params.addRequiredParam<UserObjectName>("eos_liq", "Liquid equation of state");
    params.addRequiredParam<UserObjectName>("eos_gas", "Gas equation of state");
    // Function computing the saturation tempature:
    //params.addParam<FunctionName>("saturation_temperature", "Function computing the saturation temperature.");
    return params;
}

InterfacialRelaxationTransfer::InterfacialRelaxationTransfer(const std::string & name, InputParameters parameters) :
    Material(name, parameters),
    // Boolean for mass and heat transfers:
    _isMassOn(getParam<bool>("isMassOn")),
    _isHeatOn(getParam<bool>("isHeatOn")),
    // Aux variable for liquid phase:
    _vel_x_l(coupledValue("velocity_x_liq")),
    _vel_y_l(_dim>=2 ? coupledValue("velocity_y_liq") : _zero),
    _vel_z_l(_dim==3 ? coupledValue("velocity_z_liq") : _zero),
    _pressure_l(coupledValue("pressure_liq")),
    _rho_l(coupledValue("density_liq")),
    _alpha_l(coupledValue("vf_liquid")),
    _grad_alpha_l(coupledGradient("vf_liquid")),
    // Aux variable for gas phase:
    _vel_x_g(coupledValue("velocity_x_gas")),
    _vel_y_g(_dim>=2 ? coupledValue("velocity_y_gas") : _zero),
    _vel_z_g(_dim==3 ? coupledValue("velocity_z_gas") : _zero),
    _pressure_g(coupledValue("pressure_gas")),
    _rho_g(coupledValue("density_gas")),
    // Declare interfacial variables:
    _Aint(declareProperty<Real>("interfacial_area")),
    _PI(declareProperty<Real>("interfacial_pressure")),
    _velI(declareProperty<RealVectorValue>("interfacial_velocity")),
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
    _Omega(declareProperty<Real>("mass_transfer")),
    // Constant:
    _Aint_max(getParam<Real>("Aint")),
    // UserObject:
    _eos_liq(getUserObject<EquationOfState>("eos_liq")),
    _eos_gas(getUserObject<EquationOfState>("eos_gas"))
    // Function computing the saturation temperature:
    //_SatTemp(isCoupled("saturation_temperature") ? getFunction("saturation_temperature") : getFunction(""))
{
}

void
InterfacialRelaxationTransfer::computeQpProperties()
{
    // Initialize velocity vectors for each phase:
    RealVectorValue _vel_l(_vel_x_l[_qp], _vel_y_l[_qp], _vel_z_l[_qp]);
    RealVectorValue _vel_g(_vel_x_g[_qp], _vel_y_g[_qp], _vel_z_g[_qp]);
    
    // Compute interfacial variables:
    if (_Aint_max == 0) {
        _PI_bar[_qp] = _pressure_l[_qp]; _velI_bar[_qp] = 0.;
        _PI[_qp] = _pressure_l[_qp]; _velI[_qp] = 0.;
        _Aint[_qp] = 0.; _P_rel[_qp] = 0.; _vel_rel[_qp] = 0.;
    }
    else {
        // Compute the speed of sound for each phase:
        Real _c2_l = _eos_liq.c2_from_p_rho(_rho_l[_qp], _pressure_l[_qp]);
        Real _c2_g = _eos_gas.c2_from_p_rho(_rho_g[_qp], _pressure_g[_qp]);
    
        // Compute the impedences for each phase:
        Real _Z_l = _rho_l[_qp] * std::sqrt(_c2_l);
        Real _Z_g = _rho_g[_qp] * std::sqrt(_c2_g);
        Real _sum_Z = _Z_l + _Z_g;
    
        // Compute unit vector based on gradient of liquid void fraction:
        Real _eps = std::sqrt(std::numeric_limits<Real>::min());
        RealVectorValue _n(_grad_alpha_l[_qp](0), _grad_alpha_l[_qp](1), _grad_alpha_l[_qp](2));
        //std::cout<<_n.size()<<std::endl;
        if (_n.size() <= 1e-8) {
            _n(0) = 0.; _n(1) = 0.; _n(2) = 0.;
        }
        else {
            _n = _n / (_n.size() + _eps); }
        //std::cout<<_n<<std::endl;
        //std::cout<<_n.size()<<std::endl;
        // Compute the average interfacial Relaxation parameters:
        _PI_bar[_qp] = ( _Z_g*_pressure_l[_qp] + _Z_l*_pressure_g[_qp] ) / _sum_Z;
        _velI_bar[_qp] = ( _Z_l*_vel_l + _Z_g*_vel_g ) / _sum_Z;
        //std::cout<<"PI_bar="<<_PI_bar[_qp]<<std::endl;
        //std::cout<<"velI_bar="<<_velI_bar[_qp]<<std::endl;
    
        // Compute interfacial Relaxation parameters:
        _PI[_qp] = _PI_bar[_qp] + _Z_l*_Z_g/_sum_Z * _n*(_vel_g-_vel_l);
        _velI[_qp] = _velI_bar[_qp] + _n * (_pressure_g[_qp]-_pressure_l[_qp])/_sum_Z;
    
        // Compute the relaxation parameters:
        _Aint[_qp] = _Aint_max*(6.75*(1-_alpha_l[_qp])*(1-_alpha_l[_qp])*_alpha_l[_qp]);
        _P_rel[_qp] = _Aint[_qp] / _sum_Z; /*(mu)*/
        _vel_rel[_qp] = 0.5*_P_rel[_qp]*_Z_g*_Z_l; /*(lambda)*/
        //std::cout<<"PI="<<_PI[_qp]<<std::endl;
        //std::cout<<"velI="<<_velI[_qp]<<std::endl;
    }
    
    // Newton solve for computing TI from PI:
    Real _temp = _eos_liq.temperature_from_p_rho(_pressure_l[_qp], _rho_l[_qp]);
    if (_isMassOn == true) {
        // Define some parameters used in the local Newton solve: (DEM paper)
        Real _A = (_eos_liq.Cp() - _eos_gas.Cp() + _eos_gas.qcoeff_prime() - _eos_liq.qcoeff_prime()) / (_eos_gas.Cp() - _eos_gas.Cv());
        Real _B = (_eos_liq.qcoeff() - _eos_gas.qcoeff()) / (_eos_gas.Cp() - _eos_gas.Cv());
        Real _C = (_eos_gas.Cp() - _eos_liq.Cp()) / (_eos_gas.Cp() - _eos_gas.Cv());
        Real _D = (_eos_liq.Cp() - _eos_liq.Cv()) / (_eos_gas.Cp() - _eos_gas.Cv());
        /*std::cout<<"A="<<_A<<std::endl;
         std::cout<<"B="<<_B<<std::endl;
         std::cout<<"C="<<_C<<std::endl;
         std::cout<<"D="<<_D<<std::endl;*/
    
        // Compute the constant residual _R:
        Real _p_term_liq = std::log(_PI[_qp]+_eos_liq.Pinf());
        Real _p_term_gas = std::log(_PI[_qp]+_eos_gas.Pinf());
        Real _R = _A + _D * _p_term_liq - _p_term_gas;
        //std::cout<<"R="<<_R<<std::endl;
        //std::cout<<"Cp_liq="<<_eos_liq.Cp()<<std::endl;
        //std::cout<<"Cp_gas="<<_eos_gas.Cp()<<std::endl;
    
        // Initialyze some values and pick a guess for the temperature:
        Real _f_norm = 1;
        Real _f = 0.0;
        Real _f_prime = 0.0;
    
        // Newton solve:
        while ( std::fabs(_f_norm) > 1e-5)
        {
            _f = _R + _B / _temp + _C * std::log(_temp);
            _f_prime = _C / _temp - _B / (_temp*_temp);
            //std::cout<<"f="<<_f<<std::endl;
            //std::cout<<"f_prime="<<_f_prime<<std::endl;
            _temp = _temp - _f / _f_prime;
            //std::cout<<"T="<<_temp<<std::endl;
            _f_norm = _f / _f_prime;
        }
    }
    // Compute the interfacial temperature and density:
    _tempI[_qp] = _temp;//_SatTemp.value(_PI[_qp], _node_coord);
    _rhoI[_qp] = _eos_liq.rho_from_p_T(_pressure_l[_qp], _tempI[_qp]);
    //std::cout<<"rhoI="<<_rhoI[_qp]<<std::endl;
    //std::cout<<"tempI="<<_tempI[_qp]<<std::endl;
    
    // Compute the heat of vaporization:
    Real _Lv = (_eos_gas.Cp()-_eos_liq.Cp())*_tempI[_qp] + (_eos_gas.qcoeff()-_eos_liq.qcoeff());
    //std::cout<<"Lv="<<_Lv<<std::endl;
    // Compute the interfacial energy for liquid and gas phases:
    _EI_liq[_qp] = _eos_liq.Cp()*_tempI[_qp] + _eos_liq.qcoeff() + 0.5*_velI[_qp]*_velI[_qp];
    _EI_gas[_qp] = _eos_gas.Cp()*_tempI[_qp] + _eos_gas.qcoeff() + 0.5*_velI[_qp]*_velI[_qp];
    //std::cout<<"EI_gas="<<_EI_gas[_qp]<<std::endl;
    //std::cout<<"EI_liq="<<_EI_liq[_qp]<<std::endl;
    // Compute the heat transfer coefficient for liquid phase:
    Real _radius = 1.;
    if (_Aint_max != 0)
        _radius = 3*(1-_alpha_l[_qp])/_Aint[_qp];
    _ht_liq[_qp] = (double)_isHeatOn*5*_eos_liq.k()/_radius;
    
    // Compute the heat transfer coefficient for gas phase:
    Real _L = 2*_radius;
    Real _Re = _rho_g[_qp]*_L*std::fabs(_vel_l.size()-_vel_g.size())/_eos_gas.visc();
    Real _Pr = _eos_gas.visc()*_eos_gas.Cp()/_eos_gas.k();
    Real _Nu = 2 + 0.6*std::pow(_Re,0.5) * std::pow(_Pr, 0.33);
    _ht_gas[_qp] = (double)_isHeatOn*_eos_gas.k()*_Nu/_L;
    //std::cout<<"ht_liq="<<_ht_liq[_qp]<<std::endl;
    //std::cout<<"ht_gas="<<_ht_gas[_qp]<<std::endl;
    
    // Compute the mass transfer coefficient between phases:
    Real _temp_liq = _eos_liq.temperature_from_p_rho(_pressure_l[_qp], _rho_l[_qp]);
    Real _temp_gas = _eos_gas.temperature_from_p_rho(_pressure_g[_qp], _rho_g[_qp]);
    _Omega[_qp] = (double)_isMassOn*( _ht_liq[_qp]*(_temp_liq-_tempI[_qp])+_ht_gas[_qp]*(_temp_gas-_tempI[_qp]) ) / _Lv;
    //std::cout<<"Omega="<<_Omega[_qp]<<std::endl;
    
}