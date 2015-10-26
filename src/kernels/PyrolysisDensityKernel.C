#include "PyrolysisDensityKernel.h"

template<>
InputParameters validParams<PyrolysisDensityKernel>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription(" rho changes with T");
  params.addRequiredCoupledVar("temperature","Coupled Variable");
  return params;
}

PyrolysisDensityKernel::PyrolysisDensityKernel(const InputParameters & parameters) :
     Kernel(parameters),
     _rhov(1448),
     _rhoc(1185),
     _B(11000),
     _ER(10000),
     _nn(2),
	 _T(coupledValue("temperature")),
    _T_num(coupled("temperature"))

{
}



 Real PyrolysisDensityKernel::computeQpResidual()
{
  return  _u_dot[_qp]*_test[_i][_qp]+Source(_u[_qp],_T[_qp]) *_test[_i][_qp];
}

Real PyrolysisDensityKernel::computeQpJacobian()
{
	Real epsi = 1E-08;
	Real source = Source(_u[_qp],_T[_qp]);
	Real source_new =  Source(_u[_qp]+epsi,_T[_qp]);
  return  _test[_i][_qp]*_phi[_j][_qp]*_du_dot_du[_qp]+(source_new-source)/epsi*_phi[_j][_qp]*_test[_i][_qp];
}

 Real  PyrolysisDensityKernel::computeQpOffDiagJacobian(unsigned int jvar)
 {


	 if ( jvar == _T_num)
	    {
		 Real epsi = 1E-08;
	     Real source = Source(_u[_qp],_T[_qp]);
	 	 Real source_new =  Source(_u[_qp],_T[_qp]+epsi);
	    return (source_new-source)/epsi*_phi[_j][_qp]*_test[_i][_qp];
	    }
	 else
	 {
		 return (0);
	 }
 }
Real PyrolysisDensityKernel::Source(Real rho,Real T)
{
	if (T>600)
	{
		if ((rho-_rhoc)<1)
		{
			return 0;
		}
		else
		{
			Real tmp = (rho-_rhoc)/_rhov;
			return  _B*exp(-_ER / T)*_rhov*pow(tmp,_nn);
		}
	}
	else
	{
		return 0;
	}
}





