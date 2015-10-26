#include "PyrolysisGasPressure.h"

template<>
InputParameters validParams<PyrolysisGasPressure>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("gas pressure change with temperature");
  params.addRequiredCoupledVar("temperature","Coupled Variable");
  params.addRequiredCoupledVar("rhog","Coupled Variable");
  return params;
}

PyrolysisGasPressure::PyrolysisGasPressure(const InputParameters & parameters) :
     Kernel(parameters),
	 _R(1),
	 _mg(1),
	 _rhog(coupledValue("rhog")),
	 _rhog_num(coupled("temperature")),
     _T(coupledValue("temperature")),
     _T_num(coupled("temperature"))
{
}

Real PyrolysisGasPressure::computeQpResidual()
{
      Real r(0);
      if (_rhog[_qp]>=0)
	  {
    	  r=_u[_qp]*_test[_i][_qp]-_rhog[_qp]*_T[_qp]*_test[_i][_qp];
          return r;
	  }
      else
      {
    	  r=_u[_qp]*_test[_i][_qp];
    	  return r;
     }
}



Real PyrolysisGasPressure::computeQpJacobian()
{
	  return (_phi[_j][_qp]*_test[_i][_qp] ) ;
}
Real  PyrolysisGasPressure::computeQpOffDiagJacobian(unsigned int jvar)
 {


//	 if ( jvar == _T_num)
//	    {
//	    return -(_phi[_j][_qp]*_test[_i][_qp]);
//	    }
//	 else
//	 {
		 return (0);
//	 }

 }




