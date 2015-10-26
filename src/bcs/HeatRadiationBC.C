
#include "HeatRadiationBC.h"
#include "libmesh/meshfree_interpolation.h"
#include "libmesh/system.h"
#include "libmesh/radial_basis_interpolation.h"

template<>
InputParameters validParams<HeatRadiationBC>()
{
  InputParameters params = validParams<IntegratedBC>();
  return params;
}

HeatRadiationBC::HeatRadiationBC(const InputParameters & parameters) :
  IntegratedBC(parameters),
  _flux(getMaterialProperty<Real>("heat_flux")),
  _flux_jacobi(getMaterialProperty<Real>("heat_flux_jacobi"))
{
}

Real HeatRadiationBC::computeQpResidual()
{
  Real flux = _flux[_qp];
  return -_test[_i][_qp]*flux;
}

Real HeatRadiationBC::computeQpJacobian()
{

  Real jacobi = _flux_jacobi[_qp];
  return -_test[_i][_qp]*_phi[_j][_qp]*jacobi;
}


