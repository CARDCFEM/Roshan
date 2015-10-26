#pragma once

#include "Kernel.h"

class PyrolysisGasSource : public Kernel
{
public:

	PyrolysisGasSource(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);
private:
    VariableValue & _rho_dt;
    VariableValue & _rhoDot_dRho;
    int _rho_num;
};

template<>
InputParameters validParams<PyrolysisGasSource>();
