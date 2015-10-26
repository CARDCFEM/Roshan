
#include "GetEpsilon.h"

template<>
InputParameters validParams<GetEpsilon>()
{
  InputParameters params = validParams<AuxKernel>();
  return params;
}

GetEpsilon::GetEpsilon(const InputParameters &parameters) :
    AuxKernel(parameters),
	_epsilon_real(getMaterialProperty<Real>("epsilon_material"))
{
}

Real GetEpsilon::computeValue()
{
	Real epsilon = _epsilon_real[_qp];
   return epsilon;
}
