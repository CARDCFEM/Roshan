#pragma once

#include "Kernel.h"

class PyrolysisSource : public Kernel
{
public:

	PyrolysisSource(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);
private:
    Real _deltaH;
    VariableValue & _rho_dt;
    VariableValue & _rhoDot_dRho;
    int _rho_num;
};

template<>
InputParameters validParams<PyrolysisSource>();
