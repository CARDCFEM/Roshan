
#include "ComputeTemperatureBar.h"
#include "FEProblem.h"
#include "libmesh/elem.h"
#include "libmesh/fe.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/plane.h"

#include "RayLine.h"
#include "MooseRandom.h"

using namespace std;

template<>
InputParameters validParams<ComputeTemperatureBar>()
{
	InputParameters params = validParams<SideUserObject>();
	params += validParams<RandomInterface>();
	params.addCoupledVar("temperature", "温度场");

	return params;
}

ComputeTemperatureBar::ComputeTemperatureBar(const InputParameters & parameters) :
	SideUserObject(parameters),
	RandomInterface(parameters, *parameters.get<FEProblem *>("_fe_problem"), parameters.get<THREAD_ID>("_tid"), false),
	_temperature(coupledValue("temperature"))
{
}

Real ComputeTemperatureBar::getTemperature_Pow4_Bar(Elem * elem)  const
{
//	cout << _temperature_pow4_bar.size() <<endl;
//	if(_temperature_pow4_bar.size() == 0) return 8.1e+09;
	if(_t_step == 1) return  8.1e+09;
	for(map<const Elem *, Real>::iterator it = _temperature_pow4_bar.begin(); it != _temperature_pow4_bar.end(); ++it)
	{
		Point p = it->first->centroid()-elem->centroid();
		if(p.size() < TOLERANCE)
			return it->second;
	}
	mooseError("未找到单元");
	return 0;
//	return _temperature_pow4_bar[elem];
//	cout << "pow4:" << _temperature_pow4_bar[elem] << endl;
}
Real ComputeTemperatureBar::getTemperature_Pow3_Bar(Elem * elem)  const
{
//	if(_temperature_pow3_bar.size() == 0) return 2.7e+07;
	if(_t_step == 1) return 2.7e+07;
	for(map<const Elem *, Real>::iterator it = _temperature_pow3_bar.begin(); it != _temperature_pow3_bar.end(); ++it)
	{
		Point p = it->first->centroid()-elem->centroid();
		if(p.size() < TOLERANCE)
			return it->second;
	}
	mooseError("未找到单元");
	return 0;
}

void ComputeTemperatureBar::execute()
{
	Real temp_pow4_bar(0);
	Real temp_pow3_bar(0);

	for(int _qp = 0; _qp < _q_point.size(); ++_qp)
	{
		temp_pow4_bar += (_JxW[_qp]*pow(_temperature[_qp],4) );
		temp_pow3_bar += (_JxW[_qp]*pow(_temperature[_qp],3) );
	}

	temp_pow4_bar /= _current_side_volume;
	temp_pow3_bar /= _current_side_volume;
	_temperature_pow4_bar[_current_side_elem] = temp_pow4_bar;
	_temperature_pow3_bar[_current_side_elem] = temp_pow3_bar;
//	cout << "pow4:xia" << _temperature_pow4_bar[_current_side_elem] << endl;
//	cout <<  "side_element_centre:" << _all_element[findi]->_elem->centroid() << findi << "    T_pow4_bar:" << temperature_pow4_bar[findi] << endl;

	cout << _temperature_pow3_bar.size() <<endl;
}
