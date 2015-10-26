#pragma once

#include "Kernel.h"

class PyrolysisEnergeKernel : public Kernel
{
public:

	PyrolysisEnergeKernel(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);
private:
      Real _cp;
      Real _k;
      Real _deltaH;
      VariableValue & _rho;
      VariableValue & _rho_dt;
      VariableValue & _rhoDot_dRho;
      unsigned int _rho_num;

};

template<>
InputParameters validParams< PyrolysisEnergeKernel>();
