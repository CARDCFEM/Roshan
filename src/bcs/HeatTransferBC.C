
#include "HeatTransferBC.h"

template<>
InputParameters validParams<HeatTransferBC>()
{
  InputParameters params = validParams<IntegratedBC>();
  params.addRequiredParam<Real>("h", " ");
  params.addRequiredParam<Real>("Ta", " ");
  return params;
}

HeatTransferBC::HeatTransferBC(const InputParameters & parameters) :
  IntegratedBC(parameters),
  _h(getParam<Real>("h")),
  _ta(getParam<Real>("Ta"))
{}

Real HeatTransferBC::computeQpResidual()
{
  return _test[_i][_qp]*_h*(_u[_qp]-_ta);
}

Real HeatTransferBC::computeQpJacobian()
{
  return _test[_i][_qp]*_h*_phi[_j][_qp];
}

