#include "../../include/kernels/PyroTempTimeDerivative.h"

template<>
InputParameters validParams<PyroTempTimeDerivative>()
{
  InputParameters params = validParams<TimeDerivative>();
  params.addClassDescription(" temperature change");
  params.addRequiredCoupledVar("rho","Coupled Variable");
  return params;
}

PyroTempTimeDerivative::PyroTempTimeDerivative(const InputParameters & parameters) :
  TimeDerivative(parameters),
   _cp(1256),
   _rho(coupledValue("rho")),
   _grad_rho(coupledGradient("rho")),
   _rho_num(coupled("rho"))
{
}

Real PyroTempTimeDerivative::computeQpResidual()
{
//	std::cout << _k[_qp] <<std::endl;
  Real r(0);
//   r = _Rho0 *_cp *_du_dot_du[_qp]*_test[_i][_qp];
   r = _rho[_qp]*_cp *_u_dot[_qp]*_test[_i][_qp];
  return r;
}

Real PyroTempTimeDerivative::computeQpJacobian()
{

	return  _rho[_qp] *_cp *  _test[_i][_qp]*_phi[_j][_qp]*_du_dot_du[_qp];
}





