#include "HeatConductionMaterial.h"
#include<iostream>
using std::cout;
using std::endl;

template<>
InputParameters validParams<HeatConductionMaterial>()
{
  InputParameters params = validParams<Material>();

  params.addRequiredCoupledVar("temperature", "Coupled Temperature");
  params.addParam<std::vector<Real> >("t_list", "The vector of temperature values for building the piecewise function");
  params.addParam<std::vector<Real> >("roe_list", "The vector of rho values for building the piecewise function");
  params.addParam<std::vector<Real> >("k_list", "The vector of thermal conductivity values for building the piecewise function");
  params.addParam<std::vector<Real> >("cp_list", "The vector of specific heat values for building the piecewise function");
  params.addParam<Real>("epsilon", 0, "epsilon");
  params.addParam<std::string>("property_file","The name of the property file");
  return params;
}

HeatConductionMaterial::HeatConductionMaterial(const InputParameters & parameters) :
    Material(parameters),
    _property_file(isParamValid("property_file") ? getParam<std::string>("property_file") : "" ),
    _temperature(coupledValue("temperature")),
    _k(declareProperty<Real>("thermal_conductivity")),
    _k_dT(declareProperty<Real>("thermal_conductivity_dT")),
    _cp(declareProperty<Real>("specific_heat")),
    _cp_dT( declareProperty<Real>("specific_heat_dT")),
    _rho(declareProperty<Real>("density")),
    _rho_dT( declareProperty<Real>("density_dT")),
	_epsilon_material( declareProperty<Real>("epsilon_material"))

{
		if (_property_file == "")
		{
			_roe_list= getParam<std::vector<Real> >("roe_list");
			_T_list = getParam<std::vector<Real> >("t_list");
			_cp_list= getParam<std::vector<Real> >("cp_list");
			_k_list = getParam<std::vector<Real> > ("k_list");
			_epsilon = getParam<Real> ("epsilon");
		}
		else
		{
			readFile();
		}
		_func_roe_T =  LinearInterpolation(_T_list, _roe_list);
		_func_k_T =  LinearInterpolation(_T_list, _k_list);
		_func_cp_T =  LinearInterpolation(_T_list, _cp_list);
}

void HeatConductionMaterial::readFile()
{
	    using namespace std;
		ifstream read_file(_property_file.c_str());
		if(!read_file.good())
		    mooseError("Error opening file '" +_property_file + "' from qc data.");

		    string line;
			getline(read_file, line);
			vector<Real> num;
			istringstream iss(line);
			Real f;
			while(iss >> f)
			{
			 num.push_back(f);
            }
			_tpoint = num[0];
			_epsilon = num[1];
			for(int i = 0; i < _tpoint; ++i)
			{
				int j=0;
				vector<Real> data;
				getline(read_file, line);
				istringstream iss(line);
				Real f;
				while(iss >> f)
				{
				 data.push_back(f);
				  j+=1;
				}
                _T_list.push_back(data[0]);
                _roe_list.push_back(data[1]);
                _cp_list.push_back(data[2]);
                _k_list.push_back(data[3]);
			}

			check();
}

void HeatConductionMaterial::computeProperties()
{
  Real epsi = 1E-08;
  for (unsigned int qp(0); qp < _qrule->n_points(); ++qp)
  {
    _k[qp] = _func_k_T.sample(_temperature[qp]);
    _k_dT[qp] =  (_func_k_T.sample(_temperature[qp]+epsi) - _func_k_T.sample(_temperature[qp]-epsi))/2/epsi ;
    _cp[qp] =  _func_cp_T.sample(_temperature[qp]);
    _cp_dT[qp] = (_func_cp_T.sample(_temperature[qp]+epsi) - _func_cp_T.sample(_temperature[qp]-epsi))/2/epsi ;
    _rho[qp] = _func_roe_T.sample(_temperature[qp]);
    _rho_dT[qp] = (_func_roe_T.sample(_temperature[qp]+epsi) - _func_roe_T.sample(_temperature[qp]-epsi))/2/epsi ;
    _epsilon_material[qp]= _epsilon ;
  }
}

void HeatConductionMaterial::check()
{
	if(_tpoint != _T_list.size() )
		mooseError("温度点个数不匹配");
}
