
#pragma once

#include "RayLine.h"
#include "libmesh/elem.h"
#include <map>
using std::map;

class SideElement
{

public:
	SideElement(const Elem *elem, const Point normal, Real epsilon = 0.5, Real absorptivity = 0.5, Real diffuse_reflectivity = 0.5, Real mirrors_reflectivity = 0.5);
	~SideElement(){delete _elem;};
public:
	RayLine sendRay();
	RayLine diffuseReflectRay(RayLine * rayline, Point point);
	RayLine mirrorsReflectRay(RayLine * rayline, Point point);
	Point getSideElementNormal() {return _normal;};
	Real _epsilon;
	Real _absorptivity;
	Real _diffuse_reflectivity;
	Real _mirrors_reflectivity;
//	Real RD();
//	int MaxReflectCount;
	const Elem * _elem;

	map<SideElement*, Real> _RD;

protected:
//	const Elem * _elem;
	Point _normal;
};
