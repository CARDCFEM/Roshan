
#pragma once

#include "IntegratedBC.h"

class HeatFluxBC : public IntegratedBC
{
public:
	HeatFluxBC(const InputParameters & parameters);
protected:
	virtual Real computeQpResidual();
	const Real & _value;
};

template<>
InputParameters validParams<HeatFluxBC>();
