
#pragma once

#include "Material.h"
#include "LinearInterpolation.h"

class Function;

class PyrolysisMaterial : public Material
{
public:
	 PyrolysisMaterial(const InputParameters & parameters);

protected:
  virtual void computeProperties();

private:
  void readFile();
//  void check();

  std::string _property_file;
//  int  _tpoint;
  VariableValue & _temperature;
  MaterialProperty<Real> &_rhov;
  MaterialProperty<Real> &_rhoc;
//  MaterialProperty<Real> &_pyT;
//  MaterialProperty<Real> & _charT;
  MaterialProperty<Real> & _k;
  MaterialProperty<Real> & _k_dT;
  MaterialProperty<Real> & _cp;
  MaterialProperty<Real> & _cp_dT;
  LinearInterpolation _func_k_T;
  LinearInterpolation _func_cp_T;
  std::vector<Real> _T_list;
  std::vector<Real> _k_list;
  std::vector<Real> _cp_list;
};

template<>
InputParameters validParams<PyrolysisMaterial>();
