
#include "RayLine.h"
#include "SideElement.h"
#include "libmesh/vector_value.h"

RayLine::RayLine() :
  LineSegment(Point(0), Point(1)),
  _lenght(0)
{

}

RayLine::RayLine(const Point &p0, const Point &p1, Real lenght) :
	LineSegment(p0, p1),
	_normal( (p1-p0).unit() ),
	_p0(p0),
	_p1(p1),
	_lenght(lenght)
{
}

RayLine::RayLine(const Point &p0, const Point &normal):
	LineSegment(p0, p0+normal*100000),
	_normal(normal),
	_p0(p0),
	_p1(p0+normal*100000),
	_lenght(100000)
{
}

bool RayLine::sideIntersectedByLine(const Elem* side, Point& point)
{
	bool inter(false);
	unsigned int dim = side->dim();
	if (dim == 2)
	{
		Plane plane(*side->get_node(0), *side->get_node(1), *side->get_node(2));
		inter = intersect(plane, point);
	}
	else if (dim == 1)
	{
		RayLine side_ray(*side->get_node(0), *side->get_node(1), (*side->get_node(1)-*side->get_node(0)).size());
//		std::cout << side_ray ;
		inter = intersect(side_ray, point);
	}

	else // 1D
	{
		inter = contains_point(*side->get_node(0));
		if (inter)
		point = *side->get_node(0);
	}

	if(side->contains_point(point))
		return inter;

	else
		return false;
}

bool RayLine::sideIntersectedByLine(SideElement* side, Point& point)
{
	return sideIntersectedByLine(side->_elem, point);
}

bool RayLine::intersect (const Plane & pl, Point & intersect_p) const
{
  /**
   * There are three cases in 3D for intersection of a line and a plane
   * Case 1: The line is parallel to the plane - No intersection
   *         Numerator = non-zero
   *         Denominator = zero
   *
   * Case 2: The line is within the plane - Inf intersection
   *         Numerator = zero
   *         Denominator = zero
   *
   * Case 3: The line intersects the plane at a single point
   *         Denominator = non-zero
   */

  Point pl0 = pl.get_planar_point();
  RealVectorValue N = pl.unit_normal(_p0);
  RealVectorValue I = (_p1-_p0).unit();

  Real numerator = (pl0-_p0)*N;
  Real denominator = I*N;

  // The Line is parallel to the plane
  if (std::abs(denominator) < 1.e-10)
  {
    // The Line is on the plane
    if (std::abs(numerator) < 1.e-10)
    {
      // The solution is not unique so we'll just pick an end point for now
      intersect_p = _p0;
      return true;
    }
    return false;
  }

  Real d = numerator / denominator;

  // Make sure we haven't moved off the line segment!
  if (d + libMesh::TOLERANCE < 0 || d - libMesh::TOLERANCE > (_p1-_p0).size())
    return false;

  intersect_p = d*I + _p0;

  return true;
}

bool RayLine::intersect (const RayLine & l, Point & intersect_p) const
{
  /**
   * First check for concurance:
   *
   *
   * | x1 y1 z1 1 |
   * | x2 y2 z2 1 | = (x3 - x1) * [(x2-x1) x (x4-x3)] = 0
   * | x3 y3 z3 1 |
   * | x4 y4 z4 1 |
   *
   *
   * Solve:
   *   x = _p0 + (_p1 - _p0)*s
   *   x = l.p0 + (l._p1 - l.p0)*t
   *
   *   where
   *   a = _p1 - _p0
   *   b = l._p1 - l._p0
   *   c = l._p0 - _p0
   *
   *   s = (c x b) * (a x b) / | a x b |^2
   */
  RealVectorValue a = _p1 - _p0;
  RealVectorValue b = l._p1 - l._p0;
  RealVectorValue c = l._p0 - _p0;

  RealVectorValue v = a.cross(b);

  // Check for parallel lines
  if (std::abs(v.size()) < 1.e-10 && std::abs(c.cross(a).size()) < 1.e-10)
  {
    // TODO: The lines are co-linear but more work is needed to determine and intersection point
    //       it could be the case that the line segments don't lie on the same line or overlap only
    //       a bit
    return true;
  }

  // Check that the lines are coplanar
  Real concur = c * (a.cross(b));
  if (std::abs(concur) > 1.e-10)
    return false;

  Real s = (c.cross(b) * v) / (v*v);
  Real t = (c.cross(a) * v) / (v*v);

  // if s and t are between 0 and 1 then the Line Segments intersect
  // TODO: We could handle other case of clamping to the end of Line
  //       Segements if we want to here

  if (s >= 0 && s <= 1 && t >= 0 && t <= 1)
  {
    intersect_p = _p0 + s*a;
    return true;
  }
  return false;

  /**
   * Parameteric Equation of lines
   * _p0 + t(v0) = l._p0 + u(v1)
   *
   * Case 1: Parallel Lines
   *         v0 x v1 == 0
   *
   * Case 1a: Collinear Lines
   *         v0 x v1 == 0
   *         (l._p0 - _p0) x (_p1 - _p0) == 0
   *
   * Case 2: Intersecting Lines
   *         0 <= t <= 1
   *         0 <= u <= 1
   *
   *
   * Case 1: The lines do not intersect
   *         vleft cross vright = non-zero
   *
   * Case 2: The lines are co-linear
   *         vleft cross vright = zero
   *         vleft (Denominator) = zero
   *
   * Case 3: The line intersect at a single point
   *         vleft cross vright = zero
   *         vleft (Denominator) = non-zero
  RealVectorValue v0 = _p1 - _p0;
  RealVectorValue v1 = l._p1 - l._p0;
  RealVectorValue v2 = l._p0 - _p0;

  RealVectorValue vbot = v0.cross(v1);
  RealVectorValue vtop = v2.cross(v1);

  RealVectorValue crossed = vleft.cross(vright);

  // Case 1: No intersection
  if (std::abs(vleft.cross(vright).size()) > 1.e-10)
    return false;

  // Case 2: Co-linear (just return one of the end points)
  if (std::abs(vleft.size()) < 1.e-10)
  {
    intersect_p = _p0;
    return true;
  }

  // Case 3:

  //TODO: We could detect whether the Line Segments actually overlap
  //      instead of whether the Lines are co-linear

  Real a = vright.size()/vleft.size();
  intersect_p = _p0 + a*v0;
  return true;
     */
}
