#include "GeomVector.hh"
#include "GeomFacet2d.hh"
#include "Utilities/pointDistances.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructors, destructor.
//------------------------------------------------------------------------------
// We really don't want the default constructor, but it's required to have 
// std::vectors of these.
GeomFacet2d::
GeomFacet2d():
  mVerticesPtr(0),
  mPoints(2, 0),
  mNormal(1.0, 0.0) {
  VERIFY(false);
}

GeomFacet2d::
GeomFacet2d(const std::vector<GeomFacet2d::Vector>& vertices,
            const unsigned point1,
            const unsigned point2):
  mVerticesPtr(&vertices),
  mPoints(2),
  mNormal(vertices[point2].y() - vertices[point1].y(),
          vertices[point1].x() - vertices[point2].x()) {
  mPoints[0] = point1;
  mPoints[1] = point2;
  // REQUIRE((this->point2() - this->point1()).magnitude2() > 0.0);
  REQUIRE(fuzzyEqual((this->point2() - this->point1()).unitVector().dot(mNormal), 0.0, 1.0e-6));
  REQUIRE((this->point2() - this->point1()).unitVector().cross(mNormal).z() <= 0.0);
}

GeomFacet2d::
GeomFacet2d(const GeomFacet2d& rhs):
  mVerticesPtr(rhs.mVerticesPtr),
  mPoints(rhs.mPoints),
  mNormal(rhs.mNormal) {
}

GeomFacet2d&
GeomFacet2d::
operator=(const GeomFacet2d& rhs) {
  if (this != &rhs) {
    mVerticesPtr = rhs.mVerticesPtr;
    mPoints = rhs.mPoints;
    mNormal = rhs.mNormal;
  }
  return *this;
}

GeomFacet2d::
~GeomFacet2d() {
}

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

//-------------------------------------------------------------------------------
// Recompute our normal
//-------------------------------------------------------------------------------
void
GeomFacet2d::
computeNormal() {
  CHECK(mPoints.size() == 2);
  const auto& vertices = *mVerticesPtr;
  const auto point1 = mPoints[0];
  const auto point2 = mPoints[1];
  mNormal = Vector(vertices[point2].y() - vertices[point1].y(),
                   vertices[point1].x() - vertices[point2].x());
}

}
