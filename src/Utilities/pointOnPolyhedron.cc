//------------------------------------------------------------------------------
// pointOnPolyhedron
//
// Test if a given point is on the boundary of a polyhedron.
// You should do the bounding box check before calling this method for 
// efficiency!  We do not repeat that check here since this is meant to be 
// called by code that has already done that check.
//------------------------------------------------------------------------------
#include "pointOnPolyhedron.hh"
#include "pointOnPolygon.hh"
#include "pointInPolygon.hh"
#include "lineSegmentIntersections.hh"
#include "rotationMatrix.hh"

namespace Spheral {


bool pointOnPolyhedron(const Dim<3>::Vector& p,
                       const Dim<3>::FacetedVolume& polyhedron,
                       const double tol) {
  typedef Dim<3>::Vector Vector;
  typedef Dim<3>::FacetedVolume::Facet Facet;

  const std::vector<Vector>& vertices = polyhedron.vertices();
  const std::vector<Facet>& facets = polyhedron.facets();

  // Walk the facets and check each one until we either 
  // find one we're on or we're done.
  unsigned ifacet = 0;
  bool result = false;
  while (ifacet != facets.size() and not result) {
    const Facet& facet = facets[ifacet];

    // Check if the point is in the plane of the facet.
    if (facet.compare(p, tol) == 0) 
      result = pointInPolygon(p, vertices, facet.ipoints(), facet.normal());
    ++ifacet;
  }

  return result;
}

}
