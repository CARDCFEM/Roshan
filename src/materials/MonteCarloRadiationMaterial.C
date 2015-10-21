#include "MonteCarloRadiationMaterial.h"
#include <iostream>
#include "FEProblem.h"
#include "libmesh/elem.h"
#include "libmesh/fe.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/plane.h"

#include "ComputeTemperatureBar.h"

#include "RayLine.h"
#include "MooseRandom.h"

using namespace std;

template<>
InputParameters validParams<MonteCarloRadiationMaterial>()
{
  InputParameters params = validParams<Material>();
  params += validParams<RandomInterface>();
  params.addCoupledVar("temperature", "温度场");
  params.addRequiredParam<UserObjectName>("monte_carlo", "Monte Carlo");
  return params;
}

MonteCarloRadiationMaterial::MonteCarloRadiationMaterial(const std::string & name, InputParameters parameters) :
      Material(name, parameters),
	  RandomInterface(parameters, *parameters.get<FEProblem *>("_fe_problem"), parameters.get<THREAD_ID>("_tid"), false),
	  _current_side_elem(_assembly.sideElem()),
	  _temperature(coupledValue("temperature")),
	  _uo((getUserObject<ComputeTemperatureBar>("monte_carlo"))),
	  _flux(declareProperty<Real>("heat_flux")),
      _flux_jacobi(declareProperty<Real>("heat_flux_jacobi"))
{
}

void MonteCarloRadiationMaterial::computeQpProperties()
{
//	std::cout << "flu"  << std::endl;
//	std::cout << "热流: " << _uo.getRadiationFlux(_current_side_elem) << std::endl;
//	mooseError("d");
	Real flux =  _uo.getRadiationFlux(_current_side_elem) ;
	_flux[_qp] = flux;
	_flux_jacobi[_qp] = 0;
}
