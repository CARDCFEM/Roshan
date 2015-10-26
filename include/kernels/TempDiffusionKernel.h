#pragma once

#include "Kernel.h"

class TempDiffusionKernel : public Kernel
{
public:

	 TempDiffusionKernel(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

private:
    Real _k;
};

template<>
InputParameters validParams< TempDiffusionKernel>();
