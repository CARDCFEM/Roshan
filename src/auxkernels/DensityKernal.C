
#include "DensityKernel.h"

template<>
InputParameters validParams<DensityKernel>()
{
  InputParameters params = validParams<Kernel>();
//  params.addRequiredParam<Real>("B",  "Pre-exponential factor");
//  params.addRequiredParam<Real>("ER",  "E/R");
//  params.addRequiredCoupledVar("T",  "The temperature to be bounded");
//  params.addRequiredParam<Real>("rhov",  "The virgin density to be bounded");
//  params.addRequiredCoupledVar("rhoc",  "The char density to be bounded");
//  params.addRequiredCoupledVar("pyT",  "The py-begin temperature to be bounded");
//  params.addRequiredCoupledVar("charT",  "The char-complete temperature to be bounded");
  return params;
}

DensityKernel::DensityKernel(const InputParameters &parameters) :
    Kernel(parameters),
//	_B(getParam<Real>("B")),
//	_ER(getParam<Real>("ER")),
//	_temperature(coupledValue("T")),
   _rhov(getMaterialProperty<Real>("virgindensity")),
   _rhoc(getMaterialProperty<Real>("chardensity"))
//   _rhoc(coupledValue("rhoc")),
//   _pyT(coupledValue("pyT")),
//   _charT(coupledValue("charT"))
{
}

Real DensityKernel::computeQpResidual()
{
//	_rho=_rhov;
//   return( _rho);

	std::cout << _rhov[_qp] <<","<< _rhoc[_qp]<<	std::endl;

	return 1.;
}




