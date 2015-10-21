
#pragma once

#include "SideUserObject.h"
#include "RandomInterface.h"
#include "SideElement.h"
#include <vector>
using std::vector;
using std::map;

class RayLine;

class ComputeTemperatureBar :
public SideUserObject,
public RandomInterface
{

public:
	ComputeTemperatureBar(const InputParameters & parameters);
	ComputeTemperatureBar(const std::string & name, InputParameters parameters);
	int Which_SideelementIntersectedByLine(RayLine& ray, SideElement * sideelement_i, vector<SideElement*> sideelement_vec, Point & point);
	int Find_j_of_RDij(SideElement * sideelement_i, vector<SideElement*> sideelement_vec);
	int Find_i(const Elem * elem) const;

protected :
	virtual void initialSetup();
	virtual void initialize();
	virtual void finalize();
	virtual void execute();
	virtual void threadJoin(const UserObject & uo){};

	virtual void computeRadiationFlux();
	void computeRD();

	vector<SideElement*> _all_element;

	int _max_reflect_count;
	int _particle_count;
	Real _absorptivity;
	Real _diffuse_reflectivity;
	Real _mirrors_reflectivity;

	VariableValue &_temperature;
	vector<Real> temperature_bar;
	vector<Real> flux_radiation;

public:
	Real getTemBar(int i) const {return temperature_bar[i];}
	Real getRadiationFlux(int i)  const  {return flux_radiation[i];}
	Real getRadiationFlux(const Elem * elem)  const  {return flux_radiation[Find_i(elem)];}

};

template<>
InputParameters validParams<ComputeTemperatureBar>();
