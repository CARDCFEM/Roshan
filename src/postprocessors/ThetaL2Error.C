
#include "ThetaL2Error.h"
#include "Function.h"

template<>
InputParameters validParams<ThetaL2Error>()
{
  InputParameters params = validParams<ElementL2Error>();
  return params;
}

ThetaL2Error::ThetaL2Error(const InputParameters & parameters) :
		ElementL2Error(parameters)
{
}

Real ThetaL2Error::getValue()
{
  return std::sqrt(ElementIntegralPostprocessor::getValue());
}

Real ThetaL2Error::computeQpIntegral()
{
	Real T = _u[_qp];
	Real T1 = 300;
	Real T2 = 1300;
	Real k1 = 10;
	Real k2 = 100;
	Real theta = (T-T1) + (k2-k1)/(T2-T1)/2/k1* (T-T1)* (T-T1);
	Real diff = theta - _func.value(_t, _q_point[_qp]);

//	std::cout << theta << "  " <<_func.value(_t, _q_point[_qp]) <<std::endl; ;
	return diff*diff;
}
