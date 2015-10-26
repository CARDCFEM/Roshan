#include "RoshanApp.h"

#include "../../include/auxkernels/GetEpsilon.h"
#include "PyroRhoChangeKernel.h"
#include "PyroTempTimeDerivative.h"
#include "TempDiffusionKernel.h"
#include "ComputeTemperatureBar.h"
#include "RoshanRevision.h"

#include "Moose.h"
#include "AppFactory.h"
#include "ActionFactory.h"
#include "Syntax.h"

#include "TestIC.h"
#include "strange.h"

///////////////////////////////////////////
#include "HeatConductionTimeDerivative.h"
#include "HeatConductionKernel.h"
#include "DensityKernel.h"
#include "HeatConductionMaterial.h"
#include "HeatFluxMaterial.h"
#include "PyrolysisMaterial.h"
#include "PyrolysisSource.h"
#include "PyrolysisGasPressure.h"
#include "PyrolysisGasSource.h"
#include "PyrolysisEnergeKernel.h"
#include "PyrolysisDensityKernel.h"
#include "PyrolysisGasRhochange.h"
#include "PyrolysisGasVelocity.h"
#include "NoSourcePressure.h"
#include "GasConvection.h"
#include "HeatFluxAuxKernel.h"


#include "MonteCarloRadiationMaterial.h"

#include "HeatFaceKernel.h"

//DG BC
#include "DGDirichletBC.h"
#include "IsoThermalBC.h"
#include "HeatFluxBC.h"
#include "HeatTransferBC.h"
#include "HeatRadiationBC.h"

#include "CosAuxKernel.h"
#include "RhoTimeDerivative.h"
#include "Gaspressure.h"
//#include "OutputTest.h"

#include "AnalyticSolution.h"

#include "ThetaL2Error.h"

template<>
InputParameters validParams<RoshanApp>()
{
  InputParameters params = validParams<MooseApp>();

  params.set<bool>("use_legacy_uo_initialization") = false;
  params.set<bool>("use_legacy_uo_aux_computation") = false;
  return params;
}

RoshanApp::RoshanApp(InputParameters parameters) :
    MooseApp(parameters)
{
  srand(processor_id());

  Moose::registerObjects(_factory);
  RoshanApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  RoshanApp::associateSyntax(_syntax, _action_factory);

  printHeader();
}

RoshanApp::~RoshanApp()
{
}

void
RoshanApp::registerApps()
{
#undef  registerApp
#define registerApp(name) AppFactory::instance().reg<name>(#name)
  registerApp(RoshanApp);
#undef  registerApp
#define registerApp(name) AppFactory::instance().regLegacy<name>(#name)
}

void
RoshanApp::registerObjects(Factory & factory)
{
#undef registerObject
#define registerObject(name) factory.reg<name>(stringifyName(name))

	registerInitialCondition(TestIC);
	registerInitialCondition(Strange);

	registerUserObject(ComputeTemperatureBar);
	registerMaterial(MonteCarloRadiationMaterial);


	registerKernel(HeatConductionTimeDerivative);
	registerKernel(HeatConductionKernel);
//	registerKernel(DensityKernel);
//	registerKernel(PyroRhoChangeKernel);
//	registerKernel (PyroTempTimeDerivative);
//	registerKernel (TempDiffusionKernel);
//	registerKernel (PyrolysisSource);
//	registerKernel (PyrolysisGasPressure);
//	registerKernel (PyrolysisGasSource);
//	registerKernel (PyrolysisEnergeKernel);
//	registerKernel (PyrolysisDensityKernel);
//	registerKernel (PyrolysisGasRhochange);
//	registerKernel (PyrolysisGasVelocity);
//	registerKernel (NoSourcePressure);
//	registerKernel (GasConvection);


	registerMaterial(HeatConductionMaterial);
	registerMaterial(HeatFluxMaterial);
//	registerMaterial(PyrolysisMaterial);

	registerDGKernel(HeatFaceKernel);

	registerBoundaryCondition(DGDirichletBC);
	registerBoundaryCondition(IsoThermalBC);
	registerBoundaryCondition(HeatFluxBC);
	registerBoundaryCondition(HeatTransferBC);
	registerBoundaryCondition(HeatRadiationBC);

//	registerAux(CosAuxKernel);
//	registerAux(RhoTimeDerivative);
//	registerAux(Gaspressure);
	registerAux(GetEpsilon);
	registerAux(HeatFluxAuxKernel);

	registerFunction(AnalyticSolution);

	registerPostprocessor(ThetaL2Error);

//	registerOutput(OutputTest);
#undef registerObject
#define registerObject(name) factory.regLegacy<name>(stringifyName(name))
}

void
RoshanApp::associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
}

void RoshanApp::printHeader()
{
	std::string line("*********************************\n\n");
	Moose::out << COLOR_CYAN << line << COLOR_DEFAULT;
	Moose::out << "热传导方程连续有限元解算器 Roshan \n\n";
	Moose::out << "当前版本: " <<  COLOR_MAGENTA << ROSHAN_REVISION << COLOR_DEFAULT << std::endl;
	Moose::out << COLOR_CYAN << line << COLOR_DEFAULT;
}
