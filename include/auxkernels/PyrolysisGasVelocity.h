
#include "AuxKernel.h"

class PyrolysisGasVelocity;

template<>
InputParameters validParams<PyrolysisGasVelocity>();

class PyrolysisGasVelocity : public AuxKernel
{
public:
	PyrolysisGasVelocity(const InputParameters & parameters);

    virtual ~PyrolysisGasVelocity() {}

protected:
  virtual Real computeValue();
  Real _rhog;
  RealTensorValue _kp;
  Real _viscosity;
  Real _poro;
  int _component;
  VariableGradient & _pressure_gradient;
};


