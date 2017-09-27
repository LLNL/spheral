//------------------------------------------------------------------------------
// pointOnPolygon
//
// Test if a given point is on the boundary of a polygon.
//------------------------------------------------------------------------------
#ifndef __Spheral_pointOnPolygon__
#define __Spheral_pointOnPolygon__

#include "Geometry/Dimension.hh"

namespace Spheral {
  // 2D: Check a polygon.
  bool pointOnPolygon(const Dim<2>::Vector& p,
                      const Dim<2>::FacetedVolume& polygon,
                      const double tol = 1.0e-10);

  // 2D: Polygon as a set of vertices.
  bool pointOnPolygon(const Dim<2>::Vector& p,
                      const std::vector<Dim<2>::Vector>& vertices,
                      const std::vector<unsigned>& ipoints,
                      const double tol = 1.0e-10);

  // 3D: Polygon as a set of coplanar vertices.
  bool pointOnPolygon(const Dim<3>::Vector& p,
                      const std::vector<Dim<3>::Vector>& vertices,
                      const std::vector<unsigned>& ipoints,
                      const double tol = 1.0e-10);
}
#endif

