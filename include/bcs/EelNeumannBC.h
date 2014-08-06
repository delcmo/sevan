#ifndef EELNEUMANNBC_H
#define EELNEUMANNBC_H

#include "IntegratedBC.h"
#include "EquationOfState.h"

// Forward Declarations
class EelNeumannBC;
class EquationOfState;

template<>
InputParameters validParams<EelNeumannBC>();


/**
 * The boundary condition with specified stagnation pressure and temperature
 * A void fraction boundary has also to be included to close the boundary condition for 7eqn system
 */
class EelNeumannBC : public IntegratedBC
{

public:
  EelNeumannBC(const std::string & name, InputParameters parameters);

  virtual ~EelNeumannBC(){}

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned jvar);

    enum EFlowEquationType
    {
    CONTINUITY = 1,
    XMOMENTUM = 2,
    ENERGY = 3
    };
    
    // Eqn. name to be read from input file
    std::string _eqn_name;
    // which equation (mass/momentum/energy) this BC is acting on
    MooseEnum _eqn_type;
    
    // Booleans
    bool _isLeftBC;
    bool _isLiquid;
    bool _isDensity;

    // Variables:
    VariableValue & _area;
    // Parameters:
    Real _rho;
    Real _temp;
    Real _pressure;
    Real _vel;
    Real _alpha;
    
    // Equation of state:
    const EquationOfState & _eos;
    
};

#endif // EELNEUMANNBC_H

