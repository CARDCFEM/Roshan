#pragma once

#include "TimeDerivative.h"
#include "Material.h"


class HeatConductionTimeDerivative : public TimeDerivative
{
public:
  HeatConductionTimeDerivative(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

private:
  const MaterialProperty<Real> & _cp;
  const MaterialProperty<Real> & _cp_dT;
  const MaterialProperty<Real> & _rho;
  const MaterialProperty<Real> & _rho_dT;
};

template<>
InputParameters validParams<HeatConductionTimeDerivative>();
