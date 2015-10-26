
#include "GasConvection.h"

template<>
InputParameters validParams<GasConvection>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredCoupledVar("pressure", "The pressure field.");
  params.addRequiredParam<RealTensorValue>("tensor_kp", "The Tensor of kp");
  return params;
}

GasConvection::GasConvection(const InputParameters & parameters) :
		Kernel(parameters),
	_rhog(1000),
	_kp(getParam<RealTensorValue>("tensor_kp")),
	_viscosity(1.98e-5),
	_poro(0.3),
    _pressure_gradient(coupledGradient("pressure")),
    _pressure_gradient_num(coupled("pressure"))
{
}

Real
GasConvection::computeQpResidual()
{
  return -_rhog*_kp/_viscosity/_poro*_pressure_gradient[_qp]*_grad_u[_qp]*_test[_i][_qp];
}
Real GasConvection::computeQpJacobian()
{
  return -_rhog*_kp/_viscosity/_poro*_pressure_gradient[_qp]*_grad_phi[_j][_qp]*_test[_i][_qp];
}
Real  GasConvection::computeQpOffDiagJacobian(unsigned int jvar)
 {


	 if ( jvar == _pressure_gradient_num)
	    {
	 	  return -_rhog*_kp/_viscosity/_poro*_grad_phi[_j][_qp]*_grad_u[_qp]*_test[_i][_qp];
		  }
		 else
		  {
			 return 0;
		  }


 }

