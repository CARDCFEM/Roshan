
#pragma once

#include "Material.h"
#include "libmesh/meshfree_interpolation.h"
using std::string;
using std::vector;

class Function;

class HeatFluxMaterial : public Material
{
public:
	HeatFluxMaterial(const InputParameters & parameters);

protected:
  virtual void computeProperties();

private:
  void readFile();
  void interpolate(std::vector<Real> &qc, std::vector<Real>  &ts, std::vector<Point> &pts, Real t);
  void check();
  MaterialProperty<Real> & _flux;
  MaterialProperty<Real> & _flux_jacobi;

  VariableValue & _temperature;
  VariableValue & _epsilon;

  Real _sigma;
  Real _tw0;

  string _data_file;
  string _flux_file;
  vector<vector<Real> > _src_qc;
  vector<vector<Real> > _src_ts;
  vector<vector<Point> >_src_pts;
  vector<Real> _time_step;
  vector<string> _field_name;

  vector<InverseDistanceInterpolation<LIBMESH_DIM> *  > _idis;
  int _num_pts;
  int _num_time_step;
  std::vector<Real> _scale;
  Real _fluxcoff;
};

template<>
InputParameters validParams<HeatFluxMaterial>();
