
#include "PyrolysisGasVelocity.h"

template<>
InputParameters validParams<PyrolysisGasVelocity>()
{
  InputParameters params = validParams<AuxKernel>();
  MooseEnum component("x y z");

  params.addRequiredParam<MooseEnum>("component", component, "The desired component of velocity.");
  params.addRequiredCoupledVar("pressure", "The pressure field.");
  params.addRequiredParam<RealTensorValue>("tensor_kp", "The Tensor of kp");
//  params.addRequiredParam<VariableName>("pressure", "The pressure field.");
  return params;
}

PyrolysisGasVelocity::PyrolysisGasVelocity(const InputParameters & parameters) :
    AuxKernel(parameters),
	_rhog(10),
	_kp(getParam<RealTensorValue>("tensor_kp")),
	_viscosity(1.98e-5),
	_poro(0.3),
    _component(getParam<MooseEnum>("component")),
    _pressure_gradient(coupledGradient("pressure"))
{
}

Real
PyrolysisGasVelocity::computeValue()
{
	RealVectorValue tmp = _kp/_viscosity/_poro*_pressure_gradient[_qp];
  return -tmp(_component);
//  return -_kp/_viscosity/_poro*_pressure_gradient[_qp](_component);
}
