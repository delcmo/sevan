#ifndef EELSTAGNATIONPANDTBC_H
#define EELSTAGNATIONPANDTBC_H

#include "IntegratedBC.h"
#include "EquationOfState.h"

// Forward Declarations
class EelStagnationPandTBC;
class StiffenedGasEquationOfState;

template<>
InputParameters validParams<EelStagnationPandTBC>();


/**
 * The boundary condition with specified stagnation pressure and temperature
 * A void fraction boundary has also to be included to close the boundary condition for 7eqn system
 */
class EelStagnationPandTBC : public IntegratedBC
{

public:
  EelStagnationPandTBC(const std::string & name, InputParameters parameters);

  virtual ~EelStagnationPandTBC(){}

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned jvar);

    enum EFlowEquationType
    {
    CONTINUITY = 1,
    XMOMENTUM = 2,
    YMOMENTUM = 3,
    ZMOMENTUM = 4,
    ENERGY = 5,
    VOID_FRACTION = 6
    };

    /// Eqn. name to be read from input file
    std::string _eqn_name;
    /// which equation (mass/momentum/energy) this BC is acting on
    MooseEnum _eqn_type;

    /// Coupled variables
    VariableValue & _alrhoA;
    VariableValue & _alrhouA_n;

    /// Coupled aux variables:
    VariableValue & _press_other_phase;
    VariableValue & _area;

    /// Specified stagnation variables:
    Real _p0_bc;
    Real _T0_bc;
    Real _gamma0_bc;
    Real _alpha_bc_l;

    /// Calculated rho_0, K, etc. on the boundary:
    Real _rho0_bc;
    Real _H0_bc;
    Real _K;
    Real _H_bar;

    /// Equation of state:
    const EquationOfState & _eos;

    // Boolean phase
    bool _isLiquid;
    bool _is5EquModel;
};

#endif // EELSTAGNATIONPANDTBC_H

