
#pragma once

#include "Material.h"
#include "SideUserObject.h"
#include "RandomInterface.h"
#include "SideElement.h"
#include <vector>
using std::vector;
//using std::map;

class RayLine;
class ComputeTemperatureBar;

class MonteCarloRadiationMaterial :
		public Material,
		public RandomInterface
{
public:
	MonteCarloRadiationMaterial(const InputParameters & parameters);
//	~MonteCarloRadiationMaterial(){};
	int Which_SideelementIntersectedByLine(RayLine& ray, SideElement * sideelement_i, vector<SideElement*> sideelement_vec, Point & point);
	int Find_j_of_RDij(SideElement * sideelement_i, vector<SideElement*> sideelement_vec);
	int Find_i(const Elem * elem) const;
	Real getRadiationFlux(int i)  const  {return flux_radiation[i];}
	Real getRadiationFlux(const Elem * elem)  const  {return flux_radiation[Find_i(elem)];}
	Real getRadiationFluxJacobi(const Elem * elem)  const  {return flux_radiation_jacobi[Find_i(elem)];}
	~MonteCarloRadiationMaterial();

protected:
  virtual void computeQpProperties();

  const Elem * & _current_side_elem;

  VariableValue &_temperature;
  const ComputeTemperatureBar & _uo;
  MaterialProperty<Real> & _flux;
  MaterialProperty<Real> & _flux_jacobi;

  virtual void initialSetup();
  virtual void initialize();
  virtual void finalize();

  virtual void computeRadiationFlux();
  void computeRD();

  vector<SideElement*> _all_element;

  int _max_reflect_count;
  int _particle_count;
  Real _epsilon;
  Real _absorptivity;
  Real _diffuse_reflectivity;
  Real _mirrors_reflectivity;

//  VariableValue &_temperature;
//  vector<Real> temperature_pow4_bar;
//  vector<Real> temperature_pow3_bar;
  vector<Real> flux_radiation;
  vector<Real> flux_radiation_jacobi;

};

template<>
InputParameters validParams<MonteCarloRadiationMaterial>();
