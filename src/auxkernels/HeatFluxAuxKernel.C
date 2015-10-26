
#include "HeatFluxAuxKernel.h"

template<>
InputParameters validParams<HeatFluxAuxKernel>()
{
  InputParameters params = validParams<AuxKernel>();
  return params;
}

HeatFluxAuxKernel::HeatFluxAuxKernel(const InputParameters & parameters) :
    AuxKernel(parameters),
	_flux(getMaterialProperty<Real>("heat_flux"))
{
}

Real HeatFluxAuxKernel::computeValue()
{
  return _flux[_qp];
}
