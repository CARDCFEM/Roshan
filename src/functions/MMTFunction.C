
#include "MMTFunction.h"

template<> InputParameters validParams<MMTFunction>()
{
   InputParameters params = validParams<Function>();
   return params;
}

MMTFunction::MMTFunction(const InputParameters & parameters) :
    Function(parameters)
{
}

Real MMTFunction::value(Real t, const Point &p)
{
	Real tmp = p(0)+p(1)+p(2);
    return std::sin(2*libMesh::pi*t)*tmp*tmp*tmp;
}
