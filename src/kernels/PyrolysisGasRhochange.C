
#include "PyrolysisGasRhochange.h"

template<>
InputParameters validParams<PyrolysisGasRhochange>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Pyrolysis  Source");
  params.addRequiredCoupledVar("rho","Coupled Variable");
  return params;
}

PyrolysisGasRhochange::PyrolysisGasRhochange(const InputParameters & parameters) :
     Kernel(parameters),
	 _poro(0.3),
     _rho_dt(coupledDot("rho")),
	 _rhoDot_dRho(coupledDotDu("rho")),
     _rho_num(coupled("rho"))
{
}

Real PyrolysisGasRhochange::computeQpResidual()
{
  Real r(0);
  if (_rho_dt[_qp]<0)
  {
	  r =  _poro*_u_dot[_qp]*_test[_i][_qp]+(1-_poro)*_rho_dt[_qp]*_test[_i][_qp] ;
	  return r;
  }
  else
  {
	  r = _poro*_u_dot[_qp]*_test[_i][_qp];
	  return 0;
  }
}

Real PyrolysisGasRhochange::computeQpJacobian()
{
	return  _poro*_test[_i][_qp]*_phi[_j][_qp]*_du_dot_du[_qp];
}
Real  PyrolysisGasRhochange::computeQpOffDiagJacobian(unsigned int jvar)
 {


	 if ( jvar ==   _rho_num)
	    {
		 if (_rho_dt[_qp]<0)
		  {
	 	  return( (1-_poro)*_rhoDot_dRho[_qp]*_test[_i][_qp]*_phi[_j][_qp] );
		  }
		 else
		  {
			 return 0;
		  }
		}

	 else
	 {
		 return (0);
	 }
 }





