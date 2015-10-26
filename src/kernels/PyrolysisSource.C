#include "PyrolysisSource.h"

template<>
InputParameters validParams<PyrolysisSource>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Pyrolysis  Source");
  params.addRequiredCoupledVar("rho","Coupled Variable");
  return params;
}

PyrolysisSource::PyrolysisSource(const InputParameters & parameters) :
     Kernel(parameters),
	 _deltaH(500000),
     _rho_dt(coupledDot("rho")),
	 _rhoDot_dRho(coupledDotDu("rho")),
     _rho_num(coupled("rho"))
{
}

Real PyrolysisSource::computeQpResidual()
{
//	std::cout << _k[_qp] <<std::endl;
  Real r(0);
  r =_deltaH *  _rho_dt[_qp] ;
  return r;
}

Real PyrolysisSource::computeQpJacobian()
{
	  return( _deltaH * _rhoDot_dRho[_qp]  * _test[_i][_qp] * _phi[_j][_qp] );
}
Real  PyrolysisSource::computeQpOffDiagJacobian(unsigned int jvar)
 {


	 if ( jvar ==   _rho_num)
	    {
	 	  return( _deltaH * _rhoDot_dRho[_qp]  * _test[_i][_qp] * _phi[_j][_qp] );
	    }
	 else
	 {
		 return (0);
	 }
 }




