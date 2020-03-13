#include "GeomVector.hh"
#include "GeomFacet2d.hh"
#include "Utilities/pointDistances.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Minimum distance to a point.
//------------------------------------------------------------------------------
double
GeomFacet2d::
distance(const GeomFacet2d::Vector& p) const {
  return (p - this->closestPoint(p)).magnitude();
}

//------------------------------------------------------------------------------
// Find the point on the facet closest to the given point.
//------------------------------------------------------------------------------
GeomFacet2d::Vector
GeomFacet2d::
closestPoint(const GeomFacet2d::Vector& p) const {
  return closestPointOnSegment(p, this->point1(), this->point2());
}

//------------------------------------------------------------------------------
// Compare a set of points:
//  1 => all points above.
//  0 => points both above and below (or equal).
// -1 => all points below.
//------------------------------------------------------------------------------
int
GeomFacet2d::
compare(const std::vector<GeomFacet2d::Vector>& points,
        const double tol) const {
  if (points.size() == 0) return 0;
  int result = this->compare(points[0], tol);
  int thpt;
  for (unsigned i = 1; i != points.size(); ++i) {
    thpt = this->compare(points[i], tol);
    if (thpt != result) return 0;
  }
  return result;
}

//------------------------------------------------------------------------------
// This facet does not need to be decomposed.
//------------------------------------------------------------------------------
void
GeomFacet2d::
decompose(std::vector<std::array<Vector, 2>>& subfacets) const {
  subfacets = {{point1(), point2()}};
}

}
