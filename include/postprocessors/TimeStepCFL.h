#ifndef TIMESTEPCFL_H
#define TIMESTEPCFL_H

#include "ElementPostprocessor.h"
#include "EquationOfState.h"

class TimeStepCFL;

template<>
InputParameters validParams<TimeStepCFL>();

/**
 * The inviscid time step stability limit:
 *
 * h_e \over {|\vec u| + c}
 */
class TimeStepCFL : public ElementPostprocessor
{
public:
  TimeStepCFL(const std::string & name, InputParameters parameters);
  virtual ~TimeStepCFL();

  virtual void initialize();
  virtual void execute();
  virtual Real getValue();
  virtual void threadJoin(const UserObject & uo);

protected:

  // Coupled variables
  VariableValue & _alpha_A_l;
  VariableValue & _alpha_A_rho_l;
  VariableValue & _alpha_A_rhou_l;
  VariableValue & _alpha_A_rhoE_l;
  VariableValue & _alpha_A_rho_g;
  VariableValue & _alpha_A_rhou_g;
  VariableValue & _alpha_A_rhoE_g;

  // Equation of state
  const EquationOfState & _eos_l;
  const EquationOfState & _eos_g;

  // Parameter
  Real _cfl;
  Real _value;
};


#endif // TIMESTEPCFL_H
