
#pragma once

#include "LineSegment.h"
#include "libmesh/elem.h"
#include "libmesh/plane.h"

#include <vector>
using std::vector;

class SideElement;

class RayLine : public LineSegment
{
public:
	RayLine();
	RayLine(const Point &p0, const Point &p1, Real lenght);
	RayLine(const Point &p0, const Point &normal);
	bool sideIntersectedByLine(const Elem* side, Point& point);
	bool sideIntersectedByLine(SideElement* side, Point& point);
	bool intersect (const RayLine & l, Point & intersect_p) const;
	bool intersect (const Plane & pl, Point & intersect_p) const;

	Point _normal;
	Point _p0;
	Point _p1;
//	Point & start(){return _p0;}
	Real _lenght;
	int _count;

	friend std::ostream& operator << (std::ostream& os, const RayLine& t)
	{
		os << "RayLine: " << & t << std::endl;
		os << "start: " << t._p0 << std::endl;
		os << "end: " << t._p1 << std::endl;
		os << "normal: " << t._normal <<std::endl;
		os << std::endl;
		return os;
	}
private:
//	vector<SideElement*> path;

};
