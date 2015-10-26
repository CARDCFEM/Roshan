#pragma once

#include "Kernel.h"

class PyrolysisGasRhochange : public Kernel
{
public:

	PyrolysisGasRhochange(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);
private:
    Real _poro;
    VariableValue & _rho_dt;
    VariableValue & _rhoDot_dRho;
    int _rho_num;
};

template<>
InputParameters validParams<PyrolysisGasRhochange>();
