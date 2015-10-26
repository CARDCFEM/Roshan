
#pragma once

#include "DGKernel.h"

class HeatFaceKernel : public DGKernel
{
public:
	HeatFaceKernel(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual(Moose::DGResidualType type);
  virtual Real computeQpJacobian(Moose::DGJacobianType type);

  Real _epsilon;
  Real _sigma;
};

template<>
InputParameters validParams<HeatFaceKernel>();

