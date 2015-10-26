
#pragma once

#include "NodalBC.h"

class IsoThermalBC;

class IsoThermalBC : public NodalBC
{
public:
	IsoThermalBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

  /// The value for this BC
  const Real & _value;
};

template<>
InputParameters validParams<IsoThermalBC>();
