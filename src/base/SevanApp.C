#include "SevanApp.h"
#include "Moose.h"
#include "ElkApp.h"
#include "AppFactory.h"
// Kernels
#include "EelTimeDerivative.h"
#include "EelMass.h"
#include "EelMomentum.h"
#include "EelEnergy.h"
#include "EelArtificialVisc.h"
#include "EelCMethod.h"
#include "EelVoidFraction.h"
// Auxkernels
#include "AreaAux.h"
#include "PressureAux.h"
#include "DensityAux.h"
#include "MachNumberAux.h"
#include "VelocityAux.h"
#include "TotalEnergyAux.h"
#include "InternalEnergyAux.h"
#include "TemperatureAux.h"
#include "VoidFractionAux.h"
#include "NormVectorAux.h"
#include "VariableForDissipativeTerm.h"
// Materials
#include "ComputeViscCoeff.h"
#include "InterfacialRelaxationTransfer.h"
// BCs
#include "EelStagnationPandTBC.h"
#include "EelStaticPandTBC.h"
#include "EelWallBC.h"
#include "MassInflow.h"
// ICs
#include "ConservativeVariables1DXIC.h"
#include "ConservativeVariables1DYIC.h"
#include "ConservativeVariables2DIC.h"
//Functions
#include "AreaFunction.h"
#include "SaturationTemperature.h"
// UserObjects
#include "EquationOfState.h"
#include "JumpGradientInterface.h"
#include "SmoothFunction.h"

template<>
InputParameters validParams<SevanApp>()
{
  InputParameters params = validParams<MooseApp>();
  return params;
}

SevanApp::SevanApp(const std::string & name, InputParameters parameters) :
    MooseApp(name, parameters)
{
  srand(libMesh::processor_id());
  
  Moose::registerObjects(_factory);
  ElkApp::registerObjects(_factory);
  SevanApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  ElkApp::associateSyntax(_syntax, _action_factory);
  SevanApp::associateSyntax(_syntax, _action_factory);
}

SevanApp::~SevanApp()
{
}

void
SevanApp::registerApps()
{
  registerApp(SevanApp);
}

void
SevanApp::registerObjects(Factory & factory)
{
    // Kernels
    registerKernel(EelTimeDerivative);
    registerKernel(EelMass);
    registerKernel(EelMomentum);
    registerKernel(EelEnergy);
    registerKernel(EelArtificialVisc);
    registerKernel(EelCMethod);
    registerKernel(EelVoidFraction);
    // Auxkernels
    registerAux(AreaAux);
    registerAux(PressureAux);
    registerAux(DensityAux);
    registerAux(MachNumberAux);
    registerAux(VelocityAux);
    registerAux(TotalEnergyAux);
    registerAux(InternalEnergyAux);
    registerAux(TemperatureAux);
    registerAux(VoidFractionAux);
    registerAux(NormVectorAux);
    registerAux(VariableForDissipativeTerm);
    // Materials
    registerMaterial(ComputeViscCoeff);
    registerMaterial(InterfacialRelaxationTransfer);
    // BCs
    registerBoundaryCondition(EelStagnationPandTBC);
    registerBoundaryCondition(EelStaticPandTBC);
    registerBoundaryCondition(EelWallBC);
    registerBoundaryCondition(MassInflow);
    // ICs
    registerInitialCondition(ConservativeVariables1DXIC);
    registerInitialCondition(ConservativeVariables1DYIC);
    registerInitialCondition(ConservativeVariables2DIC);
    // Functions
    registerFunction(AreaFunction);
    registerFunction(SaturationTemperature);
    //UserObjects
    registerUserObject(EquationOfState);
    registerUserObject(JumpGradientInterface);
    registerUserObject(SmoothFunction);
}

void
SevanApp::associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
}
