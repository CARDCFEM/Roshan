
#pragma once

#include "Function.h"

class AnalyticSolution : public Function
{
public:
	AnalyticSolution(const InputParameters & parameters);
    virtual Real value(Real t, const Point & p);

protected:
    Real _alpha;
    Real _length;
};

template<> InputParameters validParams<AnalyticSolution>();
