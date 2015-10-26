#include "PyrolysisGasSource.h"

template<>
InputParameters validParams<PyrolysisGasSource>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Pyrolysis  Source");
  params.addRequiredCoupledVar("rho","Coupled Variable");
  return params;
}

PyrolysisGasSource::PyrolysisGasSource(const InputParameters & parameters) :
     Kernel(parameters),
     _rho_dt(coupledDot("rho")),
	 _rhoDot_dRho(coupledDotDu("rho")),
     _rho_num(coupled("rho"))
{
}

Real PyrolysisGasSource::computeQpResidual()
{
//	std::cout << _k[_qp] <<std::endl;
  Real r(0);
  r = _rho_dt[_qp]*_test[_i][_qp] ;
  return r;
}

Real PyrolysisGasSource::computeQpJacobian()
{
	return 0;
}
Real  PyrolysisGasSource::computeQpOffDiagJacobian(unsigned int jvar)
 {


	 if ( jvar ==   _rho_num)
	    {
	 	  return(  _rhoDot_dRho[_qp] *_test[_i][_qp]*_phi[_j][_qp] );
	    }
	 else
	 {
		 return (0);
	 }
 }





