#pragma once

#include "Kernel.h"

class GasConvection : public Kernel
{
public:

	GasConvection(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

private:
    Real _rhog;
    RealTensorValue _kp;
    Real _viscosity;
    Real _poro;
    VariableGradient & _pressure_gradient;
    int _pressure_gradient_num;
};

template<>
InputParameters validParams<GasConvection>();
