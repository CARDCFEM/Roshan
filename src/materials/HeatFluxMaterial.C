#include "HeatFluxMaterial.h"
#include "libmesh/system.h"
#include "libmesh/radial_basis_interpolation.h"
#include<iostream>
using std::cout;
using std::endl;

template<>
InputParameters validParams<HeatFluxMaterial>()
{
  InputParameters params = validParams<Material>();

  params.addRequiredCoupledVar("temperature", "Coupled Temperature");
  params.addRequiredCoupledVar("epsilon", "Coupled epsilon");
  params.addRequiredParam<Real>("sigma", "sigma");
  params.addRequiredParam<Real>("tw0", "tw0");
  params.addParam<std::string>("flux_data", "flux data from file");
  params.addParam<Real>("power", 2, "distance power");
  params.addParam<int>("interp_pts", 8, "number of interpolate points");
  params.addParam<std::vector<Real> >("scale",  "scale");
  params.addParam<Real>("fluxcoff",  "fluxcoff");
  return params;
}

HeatFluxMaterial::HeatFluxMaterial(const InputParameters & parameters) :
      Material(parameters),
    _flux(declareProperty<Real>("heat_flux")),
    _flux_jacobi(declareProperty<Real>("heat_flux_jacobi")),
    _temperature(coupledValue("temperature")),
    _epsilon(coupledValue("epsilon")),
	_sigma(getParam<Real>("sigma")),
	_tw0(getParam<Real>("tw0")),
	_flux_file(isParamValid("flux_data") ? getParam<std::string>("flux_data") : ""),
	_scale(getParam<std::vector<Real> >("scale")),
	_fluxcoff(getParam<Real>("fluxcoff"))
{
	readFile();

	Real power = getParam<Real>("power");
	int n_interp_pts = getParam<int>("interp_pts");
	_field_name.push_back("qc");
	_field_name.push_back("ts");

	for(int t = 0  ; t < _num_time_step; ++t)
	{
		vector<Real> src_vals;
		for(int i = 0; i < _num_pts; ++i)
		{
			src_vals.push_back(_src_qc[t][i]);
			src_vals.push_back(_src_ts[t][i]);
		}

		_idis.push_back(new InverseDistanceInterpolation<LIBMESH_DIM>(this->comm(), n_interp_pts, power));
		_idis[t]->set_field_variables(_field_name);
		_idis[t]->get_source_points() = _src_pts[t];
		_idis[t]->get_source_vals() = src_vals;
	}
}

void HeatFluxMaterial::readFile()
{
	using namespace std;
	ifstream qc_read(_flux_file.c_str());
    if(!qc_read.good())
	     mooseError("Error opening file '" + _flux_file + "' from  data.");

    string line;
    getline(qc_read, line);
    istringstream iss(line);
    Real f;
    vector<Real> head;
    while(iss >> f)
       head.push_back(f);

    _num_pts = head[0];
    _num_time_step = head[1];
    _src_qc.resize(_num_time_step);
    _src_ts.resize(_num_time_step);
	_src_pts.resize(_num_time_step);
    for (int t = 0; t < _num_time_step; ++t)
    {
    	vector<Real> store;
    	for(int i = 0; i < _num_pts ; ++i)
    	{
    		store.resize(_num_pts);
    	   	int j = 0;
    		getline(qc_read, line);
    		istringstream transline(line);
    		Real d;
    		while (transline >> d)
    		{
    			store.push_back(d);
    			store[j] = d;
    			j+=1;
    		}

    	    Point pp(store[0]*_scale[0],store[1]*_scale[1],store[2]*_scale[2]);
    		_src_ts[t].push_back(store[3]);
    		_src_qc[t].push_back(store[4]);
    		_src_pts[t].push_back(pp);
    	}
    	_time_step.push_back(store[7]);
	}

}


void HeatFluxMaterial::interpolate(std::vector<Real> &qc, std::vector<Real>  &ts, std::vector<Point> &pts, Real t)
{
	int lower, upper;
	if(t < _time_step.front())
	{
		lower = 0;
		upper = 0;
	}
	else if(t > _time_step.back())
	{
		lower = _num_time_step - 1;
		upper = _num_time_step - 1;
	}
	else
	{
	    std::vector<double>::iterator up = std::lower_bound(_time_step.begin(), _time_step.end(), t);
	    upper = static_cast<unsigned int>(std::distance(_time_step.begin(), up));
	    if (_time_step[upper] == t)
	      lower = upper;
	    else
	      lower = upper - 1;
	}
	Real lam = 0.5;
	if(lower != upper)
	{
		lam = (_time_step[upper] - t)/(_time_step[upper] - _time_step[lower]);
	}
	vector<Real> tgt_val_lower,  tgt_val_upper;
	_idis[lower]->interpolate_field_data(_field_name, pts, tgt_val_lower);
	_idis[upper]->interpolate_field_data(_field_name, pts, tgt_val_upper);

	vector<Real> tgt_val_current;
	for(int i  = 0; i < tgt_val_lower.size(); ++i)
	{
		tgt_val_current.push_back(lam*tgt_val_lower[i] + (1-lam)*tgt_val_upper[i]);
	}

	int num_val = _field_name.size();
	for(int ip = 0; ip < pts.size(); ++ip)
	{
		qc.push_back(tgt_val_current[num_val*ip]);
		ts.push_back(tgt_val_current[num_val*ip+1]);
	}

}

void HeatFluxMaterial::computeProperties()
{
	std::vector<Real> qc, ts;
	std::vector<Point> pts;
	for(int ip = 0; ip < _q_point.size(); ++ip)
	{
		pts.push_back(_q_point[ip]);
	}
	interpolate(qc, ts, pts, _t);

	for(_qp = 0; _qp < _qrule->n_points(); ++_qp)
	{
		Real tw = _temperature[_qp];
		Real tw4 = Utility::pow<4>(tw);
		Real tw3 = Utility::pow<3>(tw);
	    _flux[_qp] = (ts[_qp] - tw)/(ts[_qp] - _tw0)*qc[_qp] - _sigma*_epsilon[_qp]*tw4;
	    _flux_jacobi[_qp] = -qc[_qp]/(ts[_qp] - _tw0) - 4*_sigma*_epsilon[_qp]*tw3;
	}
}
