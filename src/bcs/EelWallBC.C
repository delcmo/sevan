#include "EelWallBC.h"

template<>
InputParameters validParams<EelWallBC>()
{
  InputParameters params = validParams<IntegratedBC>();

    params.addRequiredParam<std::string>("equation_name", "The name of the equation this BC is acting on");
    // Coupled variables:
    params.addRequiredCoupledVar("pressure", "fluid pressure");
    params.addRequiredCoupledVar("vf_liquid", "void fraction of the liquid");
    params.addRequiredCoupledVar("area", "area");
    // Boolean
    params.addParam<bool>("isLiquid", true, "is liquid phase or not?");

  return params;
}

EelWallBC::EelWallBC(const std::string & name, InputParameters parameters) :
    IntegratedBC(name, parameters),
    // Name of the equation:
    _eqn_name(getParam<std::string>("equation_name")),
    _eqn_type("CONTINUITY, XMOMENTUM, YMOMENTUM, ENERGY, INVALID", "INVALID"),
     // Coupled variables:
    _pressure(coupledValue("pressure")),
    _vf_liquid(coupledValue("vf_liquid")),
    _area(coupledValue("area")),
    // Boolean:
    _isLiquid(getParam<bool>("isLiquid"))
{
  _eqn_type = _eqn_name;
}

Real
EelWallBC::computeQpResidual()
{
    // Compute the void fraction of the phase:
    Real _alpha = (1-(double)_isLiquid)*(1-_vf_liquid[_qp]) + (double)_isLiquid*_vf_liquid[_qp];
    // Switch statement on equation type:
    switch (_eqn_type) {
        case CONTINUITY:
            return 0.;
            break;
        case XMOMENTUM:
            return _alpha*_area[_qp]*_pressure[_qp]*_normals[_qp](0)*_test[_i][_qp];
            break;
        case YMOMENTUM:
            return _alpha*_area[_qp]*_pressure[_qp]*_normals[_qp](1)*_test[_i][_qp];
            break;
        case ENERGY:
            return 0.;
            break;
        case VOID_FRACTION:
            return 0.;
            break;
        default:
            mooseError("The equation with name: \"" << _eqn_name << "\" is not supported in the \"EelWallBC\" type of boundary condition.");
            return 0.;
            break;
    }
}

Real
EelWallBC::computeQpJacobian()
{
  // TODO
  return 0;
}

Real
EelWallBC::computeQpOffDiagJacobian(unsigned jvar)
{
  // TODO
  return 0;
}
