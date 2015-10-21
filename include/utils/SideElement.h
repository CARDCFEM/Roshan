
#pragma once

#include "RayLine.h"
#include "libmesh/elem.h"
#include <map>
using std::map;

class SideElement
{

public:
	SideElement(const Elem *elem, const Point normal, Real absorptivity = 0.5, Real diffuse_reflectivity = 0.5, Real mirrors_reflectivity = 0.5);

public:
	RayLine sendRay();
	RayLine diffuseReflectRay(RayLine * rayline, Point point);
	RayLine mirrorsReflectRay(RayLine * rayline, Point point);
	Point getSideElementNormal() {return _normal;};
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
