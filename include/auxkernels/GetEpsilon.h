
#pragma once
#include "AuxKernel.h"

class GetEpsilon : public AuxKernel
{
public:
	GetEpsilon(const InputParameters & parameters);

protected:
  virtual Real computeValue();

  const MaterialProperty<Real> & _epsilon_real;
};

template<>
InputParameters validParams<GetEpsilon>();
