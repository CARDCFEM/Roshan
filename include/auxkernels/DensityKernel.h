
#pragma once
#include "Kernel.h"

class DensityKernel : public Kernel
{
public:
	DensityKernel(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
//  Real _rho;
//  Real _B;
//  Real _ER;
//  Real _nn;
//  Real _rhov;
  const MaterialProperty<Real> & _rhov;
  const MaterialProperty<Real> & _rhoc;
//  Real _rhoc;
//  Real _pyT;
//  Real _charT;
//  VariableValue &_temperature;
};

template<>
InputParameters validParams<DensityKernel>();
