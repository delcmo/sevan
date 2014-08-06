#ifndef EELDIRICHLETBC_H
#define EELDIRICHLETBC_H

#include "NodalBC.h"
#include "EquationOfState.h"

//Forward Declarations
class EelDirichletBC;

template<>
InputParameters validParams<EelDirichletBC>();

/**
 * Implements space-dependent Dirichlet BC.
 */
class EelDirichletBC : public NodalBC
{
public:

  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  EelDirichletBC(const std::string & name, InputParameters parameters);

  virtual ~EelDirichletBC(){}

protected:

  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

    enum EFlowEquationType
    {
        VOIDFRACTION = 0,
        CONTINUITY = 1,
        XMOMENTUM = 2,
        ENERGY = 3
    };

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

#endif // EELDIRICHLETBC_H
