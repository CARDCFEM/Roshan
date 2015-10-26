#include "NoSourcePressure.h"

template<>
InputParameters validParams<NoSourcePressure>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Pressure without Source");
  return params;
}

NoSourcePressure::NoSourcePressure(const InputParameters & parameters) :
     Kernel(parameters),
	 _kp( 8.968e-9),
	 _viscosity(1.98e-5),
	 _poro(0.3)
{
}

Real NoSourcePressure::computeQpResidual()
{
//	std::cout << _k[_qp] <<std::endl;
  Real r(0);
  r = _kp/_viscosity/_poro*_grad_u[_qp]*_grad_test[_i][_qp];
  return r;
}

Real NoSourcePressure::computeQpJacobian()
{
	  return _kp/_viscosity/_poro*_grad_phi[_j][_qp] * _grad_test[_i][_qp];
}





