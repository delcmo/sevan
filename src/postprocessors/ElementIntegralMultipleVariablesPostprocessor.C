/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "ElementIntegralMultipleVariablesPostprocessor.h"

template<>
InputParameters validParams<ElementIntegralMultipleVariablesPostprocessor>()
{
  InputParameters params = validParams<ElementIntegralPostprocessor>();
    // Variable:
    params.addRequiredParam<VariableName>("variable", "The name of the variable that this object operates on");
    // Output type
    params.addRequiredParam<std::string>("output_type", "Output type: rho*c*c, rho*c*vel or rho*vel*vel.");
    // Conservative variables:
    params.addRequiredCoupledVar("alA", "alpha*A");
    params.addRequiredCoupledVar("alrhoA", "alpha*rho*A");
    params.addRequiredCoupledVar("alrhouA_x", "alpha*rho*u*A");
    params.addCoupledVar("alrhouA_y", "alpha*rho*v*A");
    params.addRequiredCoupledVar("alrhoEA", "alpha*rho*E*A");
    // Auxkernel variable:
    params.addRequiredCoupledVar("area", "area");
    // Equation of state:
    params.addRequiredParam<UserObjectName>("eos", "Equation of state");
    // Boolean for output:
    params.addParam<bool>("isLiquid", true, "boolean for phase: true-> liquid by default");
  return params;
}

ElementIntegralMultipleVariablesPostprocessor::ElementIntegralMultipleVariablesPostprocessor(const std::string & name, InputParameters parameters) :
    ElementIntegralPostprocessor(name, parameters),
    MooseVariableInterface(parameters, false),
    // Function Mach number:
    _output_name(getParam<std::string>("output_type")),
    _output_type("RHOVEL2, RHOCVEL, RHOC2, INVALID", _output_name),
    // Variable:
    _var(_subproblem.getVariable(_tid, parameters.get<VariableName>("variable"))),
    // Conservative variables
    _alA(coupledValue("alA")),
    _alrhoA(coupledValue("alrhoA")),
    _alrhouA_x(coupledValue("alrhouA_x")),
    _alrhouA_y(_mesh.dimension()>=2 ? coupledValue("alrhouA_y"): _zero),
    _alrhoEA(coupledValue("alrhoEA")),
    // Auxkernel variable:
    _area(coupledValue("area")),
    // User Objects for eos
    _eos(getUserObject<EquationOfState>("eos")),
    // Boolean
    _isLiquid(getParam<bool>("isLiquid"))
{
  addMooseVariableDependency(mooseVariable());
}

Real
ElementIntegralMultipleVariablesPostprocessor::computeQpIntegral()
{
    // Compute the phase void fraction:
    Real alpha = _isLiquid ? _alA[_qp]/_area[_qp] : std::max(1. - _alA[_qp]/_area[_qp], 0.);
//    alpha = alpha < 0 ? alpha : 0.;
    
    // Compute the pressure from the equation of state and the conservative variables:
    Real rho = _alrhoA[_qp] / (alpha * _area[_qp]);
    RealVectorValue vel(_alrhouA_x[_qp]/_alrhoA[_qp], _alrhouA_y[_qp]/_alrhoA[_qp]);
    Real rhoE = _alrhoEA[_qp] / _alA[_qp];
    Real pressure = _eos.pressure(rho, vel.size(), rhoE);
    
    // Compute the speed of sound:
    Real c2 = _eos.c2_from_p_rho(rho, pressure);
    
    // Return value:
    switch (_output_type)
    {
        case RHOVEL2:
            return rho * vel.size() * vel.size();
        case RHOCVEL:
            return rho * std::sqrt(c2) * vel.size();
        case RHOC2:
            return rho * c2;
        default:
            mooseError("The output type is not supporter by this function.");
            break;
    }
}
