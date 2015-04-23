#pragma once

#include "Kernel.h"

class RhoTestKernel : public Kernel
{
public:

 RhoTestKernel(const std::string & name, InputParameters parameters);
protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);
   Real Source(Real rho,Real T);
private:
    Real _rhov;
    Real _rhoc;
    Real _B;
    Real _ER;
    int _nn;
    VariableValue & _T;
    unsigned int _T_num;

};

template<>
InputParameters validParams<RhoTestKernel>();
