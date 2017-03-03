//------------------------------------------------------------------------------
// pointOnPolygon
//
// Test if a given point is on the boundary of a polygon.
//------------------------------------------------------------------------------
#include "pointOnPolygon.hh"
#include "lineSegmentIntersections.hh"

namespace Spheral {

bool pointOnPolygon(const Dim<2>::Vector& p,
                    const std::vector<Dim<2>::Vector>& vertices,
                    const double tol) {
  typedef Dim<2>::Vector Vector;
  const unsigned nvertices = vertices.size();
  unsigned i = 0, j;
  bool result = false;
  while (i != nvertices and not result) {
    j = (i + 1) % nvertices;
    CHECK(i < nvertices and j < nvertices);
    result = between(vertices[i], vertices[j], p, tol);
    ++i;
  }
  return result;
}

}
