#include "TestIC.h"

template<>
InputParameters validParams<TestIC>()
{
  InputParameters params = validParams<InitialCondition>();
  return params;
}

TestIC::TestIC(const InputParameters & parameters) :
    InitialCondition(parameters)
{}

Real TestIC::value(const Point & p)
{
  return std::sin(p.size());
}
