#pragma once

#include "IntegratedBC.h"
#include "libmesh/meshfree_interpolation.h"

using std::vector;
using std::string;

class HeatRadiationBC : public IntegratedBC
{
public:
	HeatRadiationBC(const InputParameters & parameters);

protected:
//  virtual void computeResidual();
//  virtual void computeJacobian();
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
//  void interpolate(std::vector<Real> &qc, std::vector<Real>  &ts, std::vector<Point> &pts, Real t);
//  void readFile();

  const MaterialProperty<Real> & _flux;
  const MaterialProperty<Real> & _flux_jacobi;
};

template<>
InputParameters validParams<HeatRadiationBC>();
