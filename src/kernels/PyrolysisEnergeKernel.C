#include "PyrolysisEnergeKernel.h"

template<>
InputParameters validParams<PyrolysisEnergeKernel>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription(" Pyrolysis Energe Equation");
  params.addRequiredCoupledVar("rho","Coupled Variable");
  return params;
}

PyrolysisEnergeKernel::PyrolysisEnergeKernel(const InputParameters & parameters) :
     Kernel(parameters),
	 _cp(1256),
     _k(0.75),
	 _deltaH(500000),
     _rho(coupledValue("rho")),
	 _rho_dt(coupledDot("rho")),
	 _rhoDot_dRho(coupledDotDu("rho")),
     _rho_num(coupled("rho"))

{
}

Real PyrolysisEnergeKernel::computeQpResidual()
{
//	std::cout << _k[_qp] <<std::endl;
  Real r(0);
  Real timedriv(0);
  Real diffusion(0);
  Real pyrosource(0);
  timedriv = (_rho[_qp]*_cp *_u_dot[_qp]*_test[_i][_qp]);
  diffusion = (_k*_grad_u[_qp]*_grad_test[_i][_qp]);
  pyrosource = (_deltaH *  _rho_dt[_qp]);
  r = timedriv + diffusion + pyrosource;
  return r;
}

Real PyrolysisEnergeKernel::computeQpJacobian()
{
	Real r(0);
	Real timedriv(0);
	Real diffusion(0);
	Real pyrosource(0);
	timedriv = (_rho[_qp] *_cp * _test[_i][_qp]*_phi[_j][_qp]*_du_dot_du[_qp]);
	diffusion = (_k*_grad_phi[_j][_qp] * _grad_test[_i][_qp]);
	pyrosource = ( _deltaH * _rhoDot_dRho[_qp]  * _test[_i][_qp] * _phi[_j][_qp] );
	r = timedriv + diffusion + pyrosource;
   return r;
}
Real  PyrolysisEnergeKernel::computeQpOffDiagJacobian(unsigned int jvar)
 {


	 if ( jvar ==   _rho_num)
	    {
	 	  return( _phi[_j][_qp]*_cp *_u_dot[_qp]*_test[_i][_qp]+_deltaH * _rhoDot_dRho[_qp]  * _test[_i][_qp] * _phi[_j][_qp] );
	    }
	 else
	 {
		 return (0);
	 }
 }




