
#include "RhoTimeDerivative.h"

template<>
InputParameters validParams<RhoTimeDerivative>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("rho","Coupled Variable");
  return params;
}

RhoTimeDerivative::RhoTimeDerivative(const InputParameters & parameters) :
    AuxKernel(parameters),
	_rho_dt(coupledDot("rho")),
    _rho_num(coupled("rho"))
{
}

Real RhoTimeDerivative::computeValue()
{
//   return abs(_rho_dt[_qp]);
   return _rho_dt[_qp];
}
