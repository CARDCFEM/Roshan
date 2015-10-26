
#pragma once
#include "AuxKernel.h"

class HeatFluxAuxKernel : public AuxKernel
{
public:
	HeatFluxAuxKernel(const InputParameters & parameters);

protected:
  virtual Real computeValue();

  const MaterialProperty<Real> & _flux;
};

template<>
InputParameters validParams<HeatFluxAuxKernel>();
