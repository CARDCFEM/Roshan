
 #include "strange.h"
#include <iostream>

template<>
InputParameters validParams<Strange>()
{
  InputParameters params = validParams<InitialCondition>();
  return params;
}

Strange::Strange(const InputParameters & parameters) :
    InitialCondition(parameters)
{}

void Strange::print()

{
  std::cout<<"initial"<<std::endl;
}
Real Strange::value(const Point & p)
{
//	print();
  return std::sin(p.size());
}




