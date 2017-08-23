//------------------------------------------------------------------------------
// pointOnPolygon
//
// Test if a given point is on the boundary of a polygon.
//------------------------------------------------------------------------------
#include "pointOnPolygon.hh"
#include "lineSegmentIntersections.hh"

using namespace std;

namespace Spheral {

bool pointOnPolygon(const Dim<2>::Vector& p,
                    const Dim<2>::FacetedVolume& polygon,
                    const double tol) {
  typedef Dim<2>::Vector Vector;
  typedef Dim<2>::FacetedVolume::Facet Facet;
  const vector<Vector>& vertices = polygon.vertices();
  const vector<Facet>& facets = polygon.facets();
  for (const Facet& facet: facets) {
    if (between(facet.point1(), facet.point2(), p, tol)) return true;
  }
  return false;
}

}
