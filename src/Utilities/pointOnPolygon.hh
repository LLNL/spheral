//------------------------------------------------------------------------------
// pointOnPolygon
//
// Test if a given point is on the boundary of a polygon.
//------------------------------------------------------------------------------
#ifndef __Spheral_pointOnPolygon__
#define __Spheral_pointOnPolygon__

#include "Geometry/Dimension.hh"

namespace Spheral {
  bool pointOnPolygon(const Dim<2>::Vector& p,
                      const std::vector<Dim<2>::Vector>& vertices,
                      const double tol = 1.0e-10);
}
#endif

