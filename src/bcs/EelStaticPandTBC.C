#include "EelStaticPandTBC.h"

template<>
InputParameters validParams<EelStaticPandTBC>()
{
  InputParameters params = validParams<IntegratedBC>();

  params.addRequiredParam<std::string>("equation_name", "The name of the equation this BC is acting on");
    // Coupled variables:
    params.addRequiredCoupledVar("vel_x", "x component of the velocity");
    params.addCoupledVar("vel_y", "y component of the velocity");
    params.addCoupledVar("area", "area aux variable");
    params.addRequiredCoupledVar("temperature", "fluid temperature");
    params.addCoupledVar("vf_liquid", "liquid void fraction");
    // Input parameters:
    params.addRequiredParam<Real>("p_bc", "Static pressure at the boundary");
    params.addParam<Real>("T_bc", 0.0, "Liquid static temperature at the boundary");
    params.addParam<Real>("T_bc_gas", -1, "Gas static temperature at the boundary");
    params.addParam<Real>("gamma_bc", 0.0, "inflow angle for inlet BC, [-], ignored for outlet condition");
    params.addParam<Real>("alpha_bc", 1., "LIQUID void fraction at the boundary");
    // Equation of state:
    params.addRequiredParam<UserObjectName>("eos", "The name of equation of state object to use.");
    // Boolean
    params.addParam<bool>("isLiquid", true, "is liquid or not?");
  return params;
}

EelStaticPandTBC::EelStaticPandTBC(const std::string & name, InputParameters parameters) :
    IntegratedBC(name, parameters),
    // Name of the equation:
    _eqn_name(getParam<std::string>("equation_name")),
    _eqn_type("CONTINUITY, XMOMENTUM, YMOMENTUM, ENERGY, INVALID", "INVALID"),
    // Coupled variables:
    _vel_x(coupledValue("vel_x")),
    _vel_y(_mesh.dimension()>=2 ? coupledValue("vel_y") : _zero),
    _area(isCoupled("area") ? coupledValue("area") : _zero),
    _temperature(coupledValue("temperature")),
    _alpha_l(isCoupled("vf_liquid") ? coupledValue("vf_liquid") : _zero),
    // Boundary condition parameters:
    _p_bc(getParam<Real>("p_bc")),
    _T_bc(getParam<Real>("T_bc")),
    _T_bc_gas(getParam<Real>("T_bc_gas")),
    _gamma_bc(getParam<Real>("gamma_bc")),
    _alpha_bc_l(getParam<Real>("alpha_bc")),
    // Equation of state:
    _eos(getUserObject<EquationOfState>("eos")),
    // Boolean:
    _isLiquid(getParam<bool>("isLiquid"))
{
    if (_T_bc_gas < 0) {
        _T_bc_gas = _T_bc;
    }
  _eqn_type = _eqn_name;
}

Real
EelStaticPandTBC::computeQpResidual()
{
    //std::cout << "p_bc=" << _p_bc << std::endl;
    //std::cout << "T_bc=" << _T_bc << std::endl;
    //std::cout << "gamma_bc=" << _gamma_bc << std::endl;
  //
  // Compute v dot n:
  Real _v_dot_n = _vel_x[_qp] * _normals[_qp](0) + _vel_y[_qp] * _normals[_qp](1);
    if ( _v_dot_n <0 ) // Inlet
    {
        // Compute the void fraction:
        Real _alpha = (1-(double)_isLiquid)*(1-_alpha_bc_l) + (double)_isLiquid*_alpha_bc_l;
        
        //std::cout <<"##########################" << _eqn_name << std::endl;
        //std::cout << "test=" << _test[_i][_qp] << std::endl;
        //std::cout << _normals[_qp](0) << " and " << _normals[_qp](1) << std::endl;
        //std::cout << "v_dot_n=" << _v_dot_n << std::endl;
        Real _rho_bc = _eos.rho_from_p_T(_p_bc, _T_bc);
        /*std::cout << "p_bc=" << _p_bc << std::endl;
        std::cout << "T_bc=" << _T_bc << std::endl;
        std::cout << "rho_bc=" << _rho_bc << std::endl;
        std::cout << "temperature=" << _temperature[_qp] << std::endl;
        std::cout << "vel x=" << _vel_x[_qp] << std::endl;*/
        Real _e_bc = 0; Real _rhoE_bc = 0; Real _vel_y_bc = 0;
        Real _norm_vel2 = 0;
        if (_gamma_bc != 0)
            _vel_y_bc = _vel_x[_qp]*std::tan(_gamma_bc);
        //std::cout << "vel y=" << _vel_y_bc << std::endl;
        
        switch (_eqn_type) {
            case CONTINUITY:
                return _alpha*_area[_qp]*_rho_bc*(_vel_x[_qp]*_normals[_qp](0)+_vel_y_bc*_normals[_qp](1))*_test[_i][_qp];
                break;
            case XMOMENTUM:
                return _alpha*_area[_qp]*( _rho_bc*_vel_x[_qp]*(_vel_x[_qp]+_vel_y_bc) +_p_bc )*_normals[_qp](0)*_test[_i][_qp];
                break;
            case YMOMENTUM:
                return _alpha*_area[_qp]*( _rho_bc*_vel_y_bc*(_vel_x[_qp]+_vel_y_bc) +_p_bc )*_normals[_qp](1)*_test[_i][_qp];
                break;
            case ENERGY:
                _e_bc = _eos.e_from_p_rho(_p_bc, _rho_bc);
                //std::cout << "e_bc=" << _e_bc << std::endl;
                _norm_vel2 = _vel_x[_qp]*_vel_x[_qp]+_vel_y_bc*_vel_y_bc;
                _rhoE_bc = _rho_bc*(_e_bc + 0.5*_norm_vel2);
                return _alpha*_area[_qp]*( _vel_x[_qp]*(_rhoE_bc+_p_bc)*_normals[_qp](0)+_vel_y_bc*(_rhoE_bc+_p_bc)*_normals[_qp](1) )*_test[_i][_qp];
                break;
            case VOID_FRACTION:
                return 0.;
                break;
            default:
                mooseError("The equation with name: \"" << _eqn_name << "\" is not supported in the \"EelStaticPandTBC\" type of boundary condition.");
                return 0.;
                break;
        }
    
    }
    else
    {
        // Compute the void fraction:
        Real _alpha = (1-(double)_isLiquid)*(1-_alpha_l[_qp]) + (double)_isLiquid*_alpha_l[_qp];
        
        // Compute density, internal energy and total energy at the BC:
        Real _rho_bc = _eos.rho_from_p_T(_p_bc, _temperature[_qp]);
        Real _e_bc = _eos.e_from_p_rho(_p_bc, _rho_bc);
        Real _rhoE_bc = _rho_bc*(_e_bc + 0.5*_vel_x[_qp]*_vel_x[_qp]);
        /*std::cout << "p_bc=" << _p_bc << std::endl;
        std::cout << "T_bc=" << _T_bc << std::endl;
        std::cout << "rho_bc=" << _rho_bc << std::endl;
        std::cout << "temperature=" << _temperature[_qp] << std::endl;
        std::cout << "vel x=" << _vel_x[_qp] << std::endl;
        std::cout << "vel x=" << _vel_y[_qp] << std::endl;*/
        switch (_eqn_type) {
            case CONTINUITY:
                return _alpha*_area[_qp]*_rho_bc*(_vel_x[_qp]*_normals[_qp](0)+_vel_y[_qp]*_normals[_qp](1))*_test[_i][_qp];
                break;
            case XMOMENTUM:
                return _alpha*_area[_qp]*( _rho_bc*_vel_x[_qp]*(_vel_x[_qp]+_vel_y[_qp])+_p_bc )*_normals[_qp](0)*_test[_i][_qp];
                break;
            case YMOMENTUM:
                return _alpha*_area[_qp]*( _rho_bc*_vel_y[_qp]*(_vel_x[_qp]+_vel_y[_qp])+_p_bc )*_normals[_qp](1)*_test[_i][_qp];
                break;
            case ENERGY:
                return _alpha*_area[_qp]*(_vel_x[_qp]*_normals[_qp](0)+_vel_y[_qp]*_normals[_qp](1))*(_rhoE_bc + _p_bc)*_test[_i][_qp];
                break;
            case VOID_FRACTION:
                return 0.;
                break;
            default:
                mooseError("The equation with name: \"" << _eqn_name << "\" is not supported in the \"EelStaticPandTBC\" type of boundary condition.");
                return 0.;
                break;
        }
    }
}

Real
EelStaticPandTBC::computeQpJacobian()
{
  // TODO
  return 0;
}

Real
EelStaticPandTBC::computeQpOffDiagJacobian(unsigned jvar)
{
  // TODO
  return 0;
}
