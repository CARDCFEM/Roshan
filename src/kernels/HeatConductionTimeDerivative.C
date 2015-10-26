#include "HeatConductionTimeDerivative.h"

template<>
InputParameters validParams<HeatConductionTimeDerivative>()
{
  InputParameters params = validParams<TimeDerivative>();
  return params;
}

HeatConductionTimeDerivative::HeatConductionTimeDerivative(const InputParameters & parameters) :
    TimeDerivative(parameters),
    _cp(getMaterialProperty<Real>("specific_heat")),
    _cp_dT(getMaterialProperty<Real>("specific_heat_dT")),
    _rho(getMaterialProperty<Real>("density")),
	_rho_dT(getMaterialProperty<Real>("density_dT"))
{
}

Real HeatConductionTimeDerivative::computeQpResidual()
{
  return _rho[_qp]*_cp[_qp]*_u_dot[_qp]*_test[_i][_qp];
}

Real HeatConductionTimeDerivative::computeQpJacobian()
{
  Real r(0);
  r = _rho[_qp]*(_cp[_qp])*_phi[_j][_qp]*_du_dot_du[_qp]*_test[_i][_qp];
  return r;
}
