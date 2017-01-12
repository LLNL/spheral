#include <limits>

#include "GeomVector.hh"
#include "GeomFacet3d.hh"
#include "Utilities/pointDistances.hh"
#include "Utilities/pointInPolygon.hh"
#include "Utilities/lineSegmentIntersections.hh"

namespace Spheral {

using namespace std;
using std::min;
using std::max;
using std::abs;

//------------------------------------------------------------------------------
// Compare a set of points:
//  1 => all points above.
//  0 => points both above and below (or equal).
// -1 => all points below.
//------------------------------------------------------------------------------
int
GeomFacet3d::
compare(const std::vector<GeomFacet3d::Vector>& points,
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
// Compute the position.
//------------------------------------------------------------------------------
GeomFacet3d::Vector
GeomFacet3d::
position() const {
  const unsigned n = mPoints.size();
  REQUIRE(n >= 3);
  Vector result;
  unsigned i, j;
  double circum, dl;
  for (unsigned i = 0; i != n; ++i) {
    j = (i + 1) % n;
    dl = (point(i) - point(j)).magnitude();
    result += (point(i) + point(j)) * dl;
    circum += dl;
  }
  result /= 2.0*circum;
  return result;
}

//------------------------------------------------------------------------------
// Compute the area.
//------------------------------------------------------------------------------
double
GeomFacet3d::
area() const {
  Vector vecsum;
  const Vector cent = this->position();
  unsigned i, j, npts = mPoints.size();
  for (i = 0; i != npts; ++i) {
    j = (i + 1) % npts;
    vecsum += (point(i) - cent).cross(point(j) - cent);
  }
  const double result = 0.5*vecsum.magnitude();
  ENSURE(result >= 0.0);
  return result;
}

//------------------------------------------------------------------------------
// Minimum distance to a point.
//------------------------------------------------------------------------------
double
GeomFacet3d::
distance(const GeomFacet3d::Vector& p) const {
  return (p - this->closestPoint(p)).magnitude();
}

//------------------------------------------------------------------------------
// Find the point on the facet closest to the given point.
//------------------------------------------------------------------------------
GeomFacet3d::Vector
GeomFacet3d::
closestPoint(const GeomFacet3d::Vector& p) const {

  // First check if the closest point on the plane is in the facet.  If so, 
  // that's it!
  const Vector centroid = this->position();
  const Vector nhat = this->normal();
  CHECK(fuzzyEqual(nhat.magnitude2(), 1.0, 1.0e-10));
  Vector result = closestPointOnPlane(p, centroid, nhat);
  if (pointInPolygon(result, *mVerticesPtr, mPoints, nhat)) return result;

  // Nope, so look for the point around the circumference of the facet that is
  // closest.
  unsigned i, j, npts = mPoints.size();
  double r2, minr2 = numeric_limits<double>::max();
  Vector potential;
  for (i = 0; i != npts; ++i) {
    j = (i + 1) % npts;
    potential = closestPointOnSegment(p, point(i), point(j));
    r2 = (potential - p).magnitude2();
    if (r2 < minr2) {
      result = potential;
      minr2 = r2;
    }
  }

  // That's it.
  return result;
}

//------------------------------------------------------------------------------
// Check if the facet is planar.
//------------------------------------------------------------------------------
bool
GeomFacet3d::
planar() const {
  const auto nv = mPoints.size();
  if (nv == 3) return true;       // Simple if this is a triangular facet!
  const auto a = this->area();
  const auto tol = std::max(1.0e-15, 1.0e-8*a);
  const auto& verts = *mVerticesPtr;
  unsigned i = 1;
  while (i < nv - 1 and
         abs(abs((verts[mPoints[i]] - verts[mPoints[0]]).cross(verts[mPoints[(i+1)%nv]] - verts[mPoints[0]]).dot(mNormal)) - 1.0) < tol) ++i;
  return (i == nv - 1);
}

}
