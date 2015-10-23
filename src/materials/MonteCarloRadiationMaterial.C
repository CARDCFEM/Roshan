#include "MonteCarloRadiationMaterial.h"
#include <iostream>
#include "FEProblem.h"
#include "libmesh/elem.h"
#include "libmesh/fe.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/plane.h"

#include "ComputeTemperatureBar.h"

#include "RayLine.h"
#include "MooseRandom.h"

using namespace std;

template<>
InputParameters validParams<MonteCarloRadiationMaterial>()
{
  InputParameters params = validParams<Material>();
  params += validParams<RandomInterface>();
  params.addCoupledVar("temperature", "温度场");
  params.addRequiredParam<UserObjectName>("monte_carlo", "Monte Carlo");

  params.addParam<int> ("max_reflect_count", 10, "最大反射次数");
  params.addParam<int> ("particle_count", 1000000, "单元发射粒子数");
  params.addParam<Real> ("epsilon", 0.5, "发射率");
  params.addParam<Real> ("absorptivity", 0.5, "吸收率");
  params.addParam<Real> ("diffuse_reflectivity", 0.5, "漫反射百分比");
  params.addParam<Real> ("mirrors_reflectivity", 0.5, "镜反射百分比");
  return params;
}

MonteCarloRadiationMaterial::MonteCarloRadiationMaterial(const InputParameters & parameters) :
      Material(parameters),
	  RandomInterface(parameters, *parameters.get<FEProblem *>("_fe_problem"), parameters.get<THREAD_ID>("_tid"), false),
	  _current_side_elem(_assembly.sideElem()),
	  _temperature(coupledValue("temperature")),
	  _uo((getUserObject<ComputeTemperatureBar>("monte_carlo"))),
	  _flux(declareProperty<Real>("heat_flux")),
	  _flux_jacobi(declareProperty<Real>("heat_flux_jacobi")),

	  _max_reflect_count(getParam<int> ("max_reflect_count")),
	  _particle_count(getParam<int> ("particle_count")),
	  _epsilon(getParam<Real> ("epsilon")),
	  _absorptivity(getParam<Real> ("absorptivity")),
	  _diffuse_reflectivity(getParam<Real> ("diffuse_reflectivity")),
	  _mirrors_reflectivity(getParam<Real> ("mirrors_reflectivity"))
{
}

void MonteCarloRadiationMaterial::initialSetup()
{
	vector<BoundaryName> boundary = getParam<std::vector<BoundaryName> >("boundary");
	std::set<BoundaryID> boundary_ids;
	for(vector<BoundaryName>::iterator it = boundary.begin(); it != boundary.end(); ++it)
	{
		BoundaryID id = _mesh.getBoundaryID(*it);
		boundary_ids.insert(id);
	}

	MeshBase & mesh = _mesh.getMesh();
	const BoundaryInfo &bnd_info = mesh.get_boundary_info();
	MeshBase::const_element_iterator   el  = mesh.active_elements_begin();
	const MeshBase::const_element_iterator end_el = mesh.active_elements_end();
	for ( ; el != end_el ; ++el)
	{
		const Elem *elem = *el;

		for (unsigned int side=0; side < elem->n_sides(); ++side)
		{
			if (elem->neighbor(side))
				continue;

			Elem *elem_side = elem->build_side(side).release();
			int bnd_id = bnd_info.boundary_id(elem, side);
			if(find(boundary_ids.begin(), boundary_ids.end(), bnd_id) == boundary_ids.end())
				continue;

			unsigned int dim = _mesh.dimension();
			FEType fe_type(Utility::string_to_enum<Order>("CONSTANT"), Utility::string_to_enum<FEFamily>("MONOMIAL"));
			FEBase * _fe_face = (FEBase::build(dim, fe_type)).release();
			QGauss * _qface = new QGauss(dim-1, FIRST);
			_fe_face->attach_quadrature_rule(_qface);
			_fe_face->reinit(elem, side);
			const std::vector<Point> normals = _fe_face->get_normals();

			_all_element.push_back(new SideElement(elem_side, normals[0], _epsilon, _absorptivity, _diffuse_reflectivity, _mirrors_reflectivity));

			delete _fe_face;
			delete _qface;
		}
	}
	flux_radiation.resize(_all_element.size(), 0);
	flux_radiation_jacobi.resize(_all_element.size(), 0);
	computeRD();
}

void MonteCarloRadiationMaterial::initialize()
{
//	computeRadiationFlux();
}

void MonteCarloRadiationMaterial::finalize()
{
	computeRadiationFlux();
}

MonteCarloRadiationMaterial::~MonteCarloRadiationMaterial()
{
	for(int i=0; i<_all_element.size(); i++)
	{
//		if(_all_element[i])
			delete _all_element[i];
	}
}

void MonteCarloRadiationMaterial::computeRD()
{
	for(int ii  = 0; ii < _all_element.size(); ii++)
	{
		for (int i=0;i<_all_element.size();i++)
		{
			_all_element[ii]->_RD[ _all_element[i] ]=0;
		}

		SideElement * cse = _all_element[ii];

		for (int j=0;j<_particle_count;j++)
		{
			int j_of_RDij=-1;

			j_of_RDij=Find_j_of_RDij(cse, _all_element);

			if (j_of_RDij == -1)
				continue;

			else
				_all_element[ii]->_RD[ _all_element[j_of_RDij] ]=_all_element[ii]->_RD[ _all_element[j_of_RDij] ]+1.0;
		}

//		cout << endl << "单元计算结果：" << endl;
//		cout << "当前单元中心点：" << _all_element[ii]->_elem->centroid() <<endl;

		for (int i=0;i<_all_element.size();i++)
		{
			_all_element[ii]->_RD[ _all_element[i] ]=_all_element[ii]->_RD[ _all_element[i] ]/_particle_count;
//			cout << "side_element_centre:" << _all_element[i]->_elem->centroid() << "        RD:" << _all_element[ii]->_RD[ _all_element[i] ] << endl;
		}
//		mooseError("产生随机位置时不支持的网格形状：");
	}
}

int MonteCarloRadiationMaterial::Find_i(const Elem * elem) const
{
	int findi=-1;

	for (int i=0;i<_all_element.size();i++)
	{
		if( (_all_element[i]->_elem->centroid()-elem->centroid()).size()<TOLERANCE )
		{
			findi=i;
			break;
		}
	}
	return findi;
}

int MonteCarloRadiationMaterial::Which_SideelementIntersectedByLine(RayLine& ray, SideElement * sideelement_i, vector<SideElement*> sideelement_vec, Point & point)
{
	int j_max=sideelement_vec.size();
	int j_wanted=-1;
	Point pp=ray._p0;
	point=ray._p1;

	for(int j=0; j<j_max; j++)
	{
//		if( (sideelement_vec[j]->_elem->centroid()-sideelement_i->_elem->centroid()).size()<TOLERANCE )
		if( ((sideelement_vec[j]->getSideElementNormal())*(sideelement_i->getSideElementNormal()) > TOLERANCE)  )
			continue;

		else if(!(ray.sideIntersectedByLine(sideelement_vec[j]->_elem,pp)))
			continue;

		else if((pp-ray._p0).size()>(point-ray._p0).size())
			continue;

		else
		{
			j_wanted=j;
			point=pp;
		}
	}
	return j_wanted;
}

int MonteCarloRadiationMaterial::Find_j_of_RDij(SideElement * sideelement_i, vector<SideElement*> sideelement_vec)
{
	int j=0;
	int j_of_RDij=-1;
	int k=0;
	bool charge=true;
	SideElement * current_elem= sideelement_i;
	RayLine rayline_in;
	RayLine rayline_out;
	RayLine* ray_in=&rayline_in;
	RayLine* ray_out=&rayline_out;
	Point p(0,0,0);

	rayline_in=(*current_elem).sendRay();
	while (charge && (k < _max_reflect_count) )
	{
		j=Which_SideelementIntersectedByLine( rayline_in, current_elem, sideelement_vec, p);

		if(j==-1)
			return -1;

		else if(MooseRandom::rand()<sideelement_vec[j]->_absorptivity)
		{
//			cout << "Absorptivity" << endl;
			charge=false;
			j_of_RDij=j;
			break;
		}

		else if(MooseRandom::rand()<sideelement_vec[j]->_diffuse_reflectivity)
		{
//			cout << "Diffuse_Reflectivity" << endl;
			rayline_out=sideelement_vec[j]->diffuseReflectRay(ray_in,p);
			current_elem=sideelement_vec[j];
			rayline_in=rayline_out;
			j_of_RDij=j;
			k++;
			continue;
		}

		else
		{
//			cout << "Mirrors_ReflectRay" << endl;
			rayline_out=sideelement_vec[j]->mirrorsReflectRay(ray_in,p);
			current_elem=sideelement_vec[j];
			rayline_in=rayline_out;
			j_of_RDij=j;
			k++;
			continue;
		}
	}

	if(!charge)
		return j_of_RDij;

	else if(k == _max_reflect_count)
		return j_of_RDij;

	else
		return -1;
}

void MonteCarloRadiationMaterial::computeRadiationFlux()
{
	Real epsilon=5.67e-8;
	Real Flux_Rad=0.0;

	for (int i=0;i<_all_element.size();i++)
	{
		Flux_Rad=0.0;
		for (int j=0;j<_all_element.size();j++)
		{
//			Flux_Rad += (_all_element[j]->_RD[ _all_element[i] ])*_all_element[j]->_elem->volume()*_absorptivity*temperature_pow4_bar[j];
			Flux_Rad += (_all_element[j]->_RD[ _all_element[i] ])*_all_element[j]->_elem->volume()*_absorptivity*_uo.getTemperature_Pow4_Bar(_all_element[j]->_elem);
		}

//		flux_radiation[i]= epsilon*Flux_Rad/_all_element[i]->_elem->volume()-epsilon*_absorptivity*temperature_pow4_bar[i];
		flux_radiation[i]= epsilon*Flux_Rad/_all_element[i]->_elem->volume();
//		flux_radiation_jacobi[i]= 4*epsilon*(_all_element[i]->_RD[ _all_element[i] ])*_absorptivity*temperature_pow3_bar[i];
		flux_radiation_jacobi[i]= 4*epsilon*(_all_element[i]->_RD[ _all_element[i] ])*_absorptivity*_uo.getTemperature_Pow3_Bar(_all_element[i]->_elem);
//		cout << "side_element_centre:" << _all_element[i]->_elem->centroid() << i << "      Flux:" << flux_radiation[i]  << endl;
	}
}

void MonteCarloRadiationMaterial::computeQpProperties()
{
//	std::cout << "flu"  << std::endl;
//	std::cout << "热流: " << _uo.getRadiationFlux(_current_side_elem) << std::endl;
//	mooseError("d");
	Real sigma=5.67e-8;
	Real epsilon=1.0;

	int i=Find_i(_current_side_elem);

//	Real flux =  _uo.getRadiationFlux(_current_side_elem) ;
//	Real flux_jacobi =  _uo.getRadiationFluxJacobi(_current_side_elem) ;
	Real flux = flux_radiation[i];
	Real flux_jacobi =  flux_radiation_jacobi[i];

	_flux[_qp] = flux-sigma*epsilon*pow(_temperature[_qp],4);
	_flux_jacobi[_qp] = flux_jacobi-4*sigma*epsilon*pow(_temperature[_qp], 3);
//	cout << "side_element_centre:" << _current_side_elem->centroid() << "      Flux:" << _flux[_qp]  << endl;
}
