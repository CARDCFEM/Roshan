
#include "DGDirichletBC.h"
#include "Function.h"

#include "libmesh/numeric_vector.h"

#include <cmath>

template<>
InputParameters validParams<DGDirichletBC>()
{
  InputParameters params = validParams<IntegratedBC>();
  params.addParam<Real>("value", 0.0, "The value the variable should have on the boundary");
  params.addParam<FunctionName>("function", "The forcing function.");
  params.addParam<Real>("sigma", 6, "Sigma");
  params.addParam<Real>("epsilon", -1, "Epsilon");

  return params;
}

DGDirichletBC::DGDirichletBC(const InputParameters & parameters) :
    IntegratedBC(parameters),
    _func(isParamValid("function") ? &getFunction("function") : NULL),
//    _func(getFunction("function")),
	_value(getParam<Real>("value")),
    _epsilon(getParam<Real>("epsilon")),
    _sigma(getParam<Real>("sigma"))
{
}

Real DGDirichletBC::computeQpResidual()
{
  const unsigned int elem_b_order = static_cast<unsigned int> (_var.getOrder());
  const double h_elem = _current_elem->volume()/_current_side_elem->volume() * 1./std::pow(elem_b_order, 2.);

  Real fn = _value;
  if(isParamValid("function"))
    fn = _func->value(_t, _q_point[_qp]);

  Real r = 0;
  r -= (_grad_u[_qp] * _normals[_qp] * _test[_i][_qp]);
  r += _epsilon * (_u[_qp] - fn) * _grad_test[_i][_qp] * _normals[_qp];
  r += _sigma/h_elem * (_u[_qp] - fn) * _test[_i][_qp];

  return r;
}

Real
DGDirichletBC::computeQpJacobian()
{
  const unsigned int elem_b_order = static_cast<unsigned int> (_var.getOrder());
  const double h_elem = _current_elem->volume()/_current_side_elem->volume() * 1./std::pow(elem_b_order, 2.);

  Real r = 0;
  r -= (_grad_phi[_j][_qp] * _normals[_qp] * _test[_i][_qp]);
  r += _epsilon * _phi[_j][_qp] * _grad_test[_i][_qp] * _normals[_qp];
  r += _sigma/h_elem * _phi[_j][_qp] * _test[_i][_qp];

  return r;
}

