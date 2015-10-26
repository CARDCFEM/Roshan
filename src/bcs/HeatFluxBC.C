#include "HeatFluxBC.h"

template<>
InputParameters validParams<HeatFluxBC>()
{
  InputParameters params = validParams<IntegratedBC>();
  params.addParam<Real>("value", 0.0, "The value of the gradient on the boundary.");
  return params;
}

HeatFluxBC::HeatFluxBC(const InputParameters & parameters) :
  IntegratedBC(parameters),
  _value(getParam<Real>("value"))
{}

Real HeatFluxBC::computeQpResidual()
{
  return -_test[_i][_qp]*_value;
}

