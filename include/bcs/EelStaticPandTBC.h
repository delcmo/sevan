#ifndef EELSTATICPANDTBC_H
#define EELSTATICPANDTBC_H

#include "IntegratedBC.h"
#include "EquationOfState.h"

// Forward Declarations
class EelStaticPandTBC;
class EquationOfState;

template<>
InputParameters validParams<EelStaticPandTBC>();

class EelStaticPandTBC : public IntegratedBC
{

public:
  EelStaticPandTBC(const std::string & name, InputParameters parameters);

  virtual ~EelStaticPandTBC(){}

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned jvar);

  enum EFlowEquationType
  {
    CONTINUITY = 0,
    XMOMENTUM = 1,
    YMOMENTUM = 2,
    ENERGY = 3,
    VOID_FRACTION = 4
  };

  // Eqn. name to be read from input file
  std::string _eqn_name;
  // which equation (mass/momentum/energy) this BC is acting on
  MooseEnum _eqn_type;
  // Coupled aux variables
  VariableValue & _vel_x;
  VariableValue & _vel_y;
  VariableValue & _area;
  VariableValue & _temperature;
    VariableValue & _alpha_l;
  // Specified pressure
  Real _p_bc;
  // Specified temperature
  Real _T_bc;
  Real _T_bc_gas;
  /// Specified inflow angle
  Real _gamma_bc;
  // Specified liquid void fraction
    Real _alpha_bc_l;
  // Equation of state
  const EquationOfState & _eos;
  // Boolean phase
    bool _isLiquid;
};

#endif // EELSTATICPANDTBC_H

