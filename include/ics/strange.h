#pragma once

#include "InitialCondition.h"

class Strange : public InitialCondition
{
public:
  Strange(const InputParameters &parameters);
  virtual void print();
  virtual Real value(const Point & p);
};

template<>
InputParameters validParams<Strange>();
