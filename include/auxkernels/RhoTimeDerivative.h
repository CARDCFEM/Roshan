#pragma once
#include "AuxKernel.h"

class RhoTimeDerivative : public AuxKernel
{
public:
	RhoTimeDerivative(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeValue();
  VariableValue & _rho_dt;
  int _rho_num;
};

template<>
InputParameters validParams<RhoTimeDerivative>();
