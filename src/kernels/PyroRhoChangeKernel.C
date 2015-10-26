#include "../../include/kernels/PyroRhoChangeKernel.h"

template<>
InputParameters validParams<PyroRhoChangeKernel>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription(" rho changes with T");
  params.addRequiredCoupledVar("temperature","Coupled Variable");
  return params;
}

 PyroRhoChangeKernel::PyroRhoChangeKernel(const InputParameters & parameters) :
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



 Real PyroRhoChangeKernel::computeQpResidual()
{
//	std::cout << _k[_qp] <<std::endl;
//  Real tmp = (_u[_qp]-_rhoc)/_rhov;
//  Real r(0);
//  r = _B*exp(-_ER/_T0)*_rhov*pow(tmp,_nn)*_test[_i][_qp];
  return  Source(_u[_qp],_T[_qp]) *_test[_i][_qp];
}

Real PyroRhoChangeKernel::computeQpJacobian()
{
	Real epsi = 1E-08;
	Real source = Source(_u[_qp],_T[_qp]);
	Real source_new =  Source(_u[_qp]+epsi,_T[_qp]);
  return (source_new-source)/epsi*_phi[_j][_qp]*_test[_i][_qp];
}

 Real  PyroRhoChangeKernel::computeQpOffDiagJacobian(unsigned int jvar)
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
Real PyroRhoChangeKernel::Source(Real rho,Real T)
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





