
#include "CosAuxKernel.h"

template<>
InputParameters validParams<CosAuxKernel>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredParam<Real>("coef",  "系数");
  params.addRequiredCoupledVar("T",  "The variable to be bounded");
  return params;
}

CosAuxKernel::CosAuxKernel(const InputParameters & parameters) :
    AuxKernel(parameters),
	_coef(getParam<Real>("coef")),
	_temperature(coupledValue("T"))
{
}

Real CosAuxKernel::computeValue()
{
	Point p = _q_point[_qp];
   return( _temperature[_qp]*_coef*cos(fabs(p(0) + p(1) + p(2))));
}
