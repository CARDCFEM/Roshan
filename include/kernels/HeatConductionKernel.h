#pragma once

#include "Kernel.h"
#include "Material.h"

class HeatConductionKernel : public Kernel
{
public:

  HeatConductionKernel(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

private:
  const MaterialProperty<Real> & _k;
  const MaterialProperty<Real> & _k_dT;
};

template<>
InputParameters validParams<HeatConductionKernel>();
