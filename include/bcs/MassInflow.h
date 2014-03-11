#ifndef MASSINFLOW_H
#define MASSINFLOW_H

#include "IntegratedBC.h"
#include "EquationOfState.h"

// Forward Declarations
class MassInflow;
class EquationOfState;

template<>
InputParameters validParams<MassInflow>();

class MassInflow : public IntegratedBC
{

public:
  MassInflow(const std::string & name, InputParameters parameters);

  virtual ~MassInflow(){}

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned jvar);

    enum EFlowEquationType
    {
      CONTINUITY = 0,
      XMOMENTUM = 1,
      YMOMENTUM = 2,
      ENERGY = 3
    };

    // Eqn. name to be read from input file
    std::string _eqn_name;
    
    // which equation (mass/momentum/energy) this BC is acting on
    MooseEnum _eqn_type;
    
    // Boolean for phase;
    bool _isLiquid;
    
    // Coupled variables:
    VariableValue & _alA;
    VariableValue & _alrhoA;
    VariableValue & _alrhouA;
    VariableValue & _alrhovA;
    VariableValue & _alrhoEA;
    
    // Coupled aux variables
    VariableValue & _area;
    
    // Specified velocity values:
    Real _momA_bc;
    Real _T_bc;
    Real _theta_bc;
    
    // Equation of state
    const EquationOfState & _eos;
};

#endif // MASSINFLOW_H

