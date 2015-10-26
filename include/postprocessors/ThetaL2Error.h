
#include "ElementL2Error.h"

class Function;

class ThetaL2Error : public ElementL2Error
{
public:
	ThetaL2Error(const InputParameters & parameters);

  virtual Real getValue();

protected:
  virtual Real computeQpIntegral();

};

template<>
InputParameters validParams<ThetaL2Error>();
