
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
	Real getTemperature_Pow4_Bar(Elem * elem)  const;
	Real getTemperature_Pow3_Bar(Elem * elem)  const;
protected :
	virtual void initialize(){};
	virtual void execute();
	virtual void threadJoin(const UserObject & uo){};
	virtual void finalize(){};
//
	VariableValue &_temperature;
	mutable map<const Elem*, Real> _temperature_pow4_bar;
	mutable map<const Elem*, Real> _temperature_pow3_bar;
};


template<>
InputParameters validParams<ComputeTemperatureBar>();
