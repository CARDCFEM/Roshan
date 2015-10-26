#pragma once

#include "IntegratedBC.h"

class DGDirichletBC : public IntegratedBC
{
public:
  DGDirichletBC(const InputParameters & parameters);
  virtual ~DGDirichletBC() {}

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

private:
  Function *_func;
  Real _value;

  Real _epsilon;
  Real _sigma;
};

template<>
InputParameters validParams<DGDirichletBC>();
