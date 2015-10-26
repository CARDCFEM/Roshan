
#include "AnalyticSolution.h"

template<> InputParameters validParams<AnalyticSolution>()
{
   InputParameters params = validParams<Function>();
   return params;
}

AnalyticSolution::AnalyticSolution(const InputParameters & parameters) :
    Function(parameters),
	_alpha(2.5E-06),
	_length(0.01)
{
}

Real AnalyticSolution::value(Real t, const Point &p)
{
	Real x = p(0);
	Real L = _length;
	Real pi = libMesh::pi;

	Real f = _alpha*t/L/L+1./3-x/L+1./2*(x/L)*(x/L);
	for(int n = 1; n < 100; n++)
	{
		f -= 2/pi/pi*(1./n/n*exp(-n*n*pi*pi*_alpha*t/L/L)*cos(n*pi*x/L));
	}

	Real k1 = 10;
	return 7.5E+05*L/k1*f;
}
