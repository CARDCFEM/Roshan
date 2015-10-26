#pragma once

#include "Kernel.h"

class PyrolysisGasPressure : public Kernel
{
public:

	PyrolysisGasPressure(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);
private:
  Real _R;
  Real _mg;
  VariableValue & _rhog;
  unsigned int _rhog_num;
  VariableValue & _T;
  unsigned int _T_num;
};

template<>
InputParameters validParams< PyrolysisGasPressure>();
