
#include "Gaspressure.h"

template<>
InputParameters validParams<Gaspressure>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addClassDescription("gas pressure change with temperature");
  params.addRequiredCoupledVar("temperature","Coupled Variable");
  params.addRequiredCoupledVar("rhog","Coupled Variable");
  return params;
}

Gaspressure::Gaspressure(const InputParameters & parameters) :
    AuxKernel(parameters),
	_R(20),
    _mg(1e-2),
    _rhog(coupledValue("rhog")),
    _rhog_num(coupled("temperature")),
	_T(coupledValue("temperature")),
    _T_num(coupled("temperature"))
{
}

Real Gaspressure::computeValue()
{
   return(_T[_qp]);
}
