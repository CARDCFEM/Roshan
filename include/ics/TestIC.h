
#pragma once

#include "InitialCondition.h"

class TestIC : public InitialCondition
{
public:
  TestIC(const InputParameters & parameters);
  virtual Real value(const Point & p);
};

template<>
InputParameters validParams<TestIC>();
