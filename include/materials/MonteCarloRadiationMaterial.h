
#pragma once

#include "Material.h"
#include "SideUserObject.h"
#include "RandomInterface.h"
#include "SideElement.h"
#include <vector>
using std::vector;
using std::map;

class RayLine;
class ComputeTemperatureBar;

class MonteCarloRadiationMaterial :
		public Material,
		public RandomInterface
{
public:
	MonteCarloRadiationMaterial(const InputParameters & parameters);
	~MonteCarloRadiationMaterial();

protected:
	int Which_SideelementIntersectedByLine(RayLine& ray, SideElement * sideelement_i, vector<SideElement*> sideelement_vec, Point & point);
	int Find_j_of_RDij(SideElement * sideelement_i, vector<SideElement*> sideelement_vec);
	int Find_i(const Elem * elem) const;
	void computeRadiationFlux();
	void computeRD();

	virtual void initialSetup();
	virtual void timestepSetup();
	virtual void computeQpProperties();
	virtual void computeProperties();

	vector<SideElement*> _all_element;
	map<Elem*, Real> _temperature_pow4_bar;
	map<Elem*, Real> _temperature_pow3_bar;
	const Elem * & _current_side_elem;
	const ComputeTemperatureBar & _uo;
	VariableValue &_temperature;
	MaterialProperty<Real> & _flux;
	MaterialProperty<Real> & _flux_jacobi;
	vector<Real> flux_radiation;
	vector<Real> flux_radiation_jacobi;


	int _max_reflect_count;
	int _particle_count;
	Real _epsilon;
	Real _absorptivity;
	Real _diffuse_reflectivity;
	Real _mirrors_reflectivity;
	Real _sigma;


};

template<>
InputParameters validParams<MonteCarloRadiationMaterial>();
