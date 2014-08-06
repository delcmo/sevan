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

#ifndef CONSERVATIVEVARIABLES1DXIC_H
#define CONSERVATIVEVARIABLES1DXIC_H

// MOOSE Includes
#include "InitialCondition.h"
#include "EquationOfState.h"
#include "AreaFunction.h"

// Forward Declarations
class ConservativeVariables1DXIC;

template<>
InputParameters validParams<ConservativeVariables1DXIC>();

/**
 * ConservativeVariables1DXIC just returns a constant value.
 */
class ConservativeVariables1DXIC : public InitialCondition
{
public:

  /**
   * Constructor: Same as the rest of the MOOSE Objects
   */
  ConservativeVariables1DXIC(const std::string & name,
            InputParameters parameters);

  /**
   * The value of the variable at a point.
   *
   * This must be overriden by derived classes.
   */
  virtual Real value(const Point & p);

private:
    // Type of ics
    bool _isDensity;
    // Function area
    Function & _area;
    // Initial conditions for left and right values:
    Real _p_left_liq;
    Real _p_right_liq;
    Real _p_left_gas;
    Real _p_right_gas;
    Real _v_left_liq;
    Real _v_right_liq;
    Real _v_left_gas;
    Real _v_right_gas;
    Real _t_left_liq;
    Real _t_right_liq;
    Real _t_left_gas;
    Real _t_right_gas;
    Real _rho_left_liq;
    Real _rho_right_liq;
    Real _rho_left_gas;
    Real _rho_right_gas;
    Real _alpha_left;
    Real _alpha_right;
    // Position of the membrane:
    Real _membrane;
    Real _length;
    // Name of the variable:
    std::string _name_var;
    // Equation of state:
    const EquationOfState & _eos;
    // Boolean for phase:
    bool _isLiquid;

};

#endif //CONSERVATIVEVARIABLES1DXIC_H
