#ifndef EELWALLBC_H
#define EELWALLBC_H

#include "IntegratedBC.h"

class EelWallBC;

template<>
InputParameters validParams<EelWallBC>();

class EelWallBC : public IntegratedBC
{

public:
  EelWallBC(const std::string & name, InputParameters parameters);

  virtual ~EelWallBC(){}

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

    /// Eqn. name to be read from input file
    std::string _eqn_name;
    /// which equation (mass/momentum/energy) this BC is acting on
    MooseEnum _eqn_type;
    // Coupled aux variables
    VariableValue & _pressure;
    VariableValue & _vf_liquid;
    VariableValue & _area;
    // Boolean for phase
    bool _isLiquid;
};

#endif // EELWALLBC_H

