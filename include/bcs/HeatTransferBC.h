#pragma once

#include "IntegratedBC.h"

class HeatTransferBC : public IntegratedBC
{
public:
  HeatTransferBC(const InputParameters & parameters);
protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

  Real _h;
  Real _ta;
};

template<>
InputParameters validParams<HeatTransferBC>();
