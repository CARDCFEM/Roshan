#include "TempDiffusionKernel.h"

template<>
InputParameters validParams<TempDiffusionKernel>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Test rho change");
  return params;
}

TempDiffusionKernel::TempDiffusionKernel(const InputParameters & parameters) :
     Kernel(parameters),
     _k(0.75)
{
}

Real TempDiffusionKernel::computeQpResidual()
{
//	std::cout << _k[_qp] <<std::endl;
  Real r(0);
  r = _k*_grad_u[_qp]*_grad_test[_i][_qp];
  return r;
}

Real TempDiffusionKernel::computeQpJacobian()
{
	  return _k*_grad_phi[_j][_qp] * _grad_test[_i][_qp];
}





