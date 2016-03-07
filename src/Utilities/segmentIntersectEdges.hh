//------------------------------------------------------------------------------
// segmentIntersectEdges.
//
// Test if a line segment (a0,a1) intersects any of the edges of a faceted
// volume.
//------------------------------------------------------------------------------
#ifndef __Spheral_segmentIntersectEdges__
#define __Spheral_segmentIntersectEdges__

#include "Geometry/Dimension.hh"

namespace Spheral {

// 2-D
bool segmentIntersectEdges(const Dim<2>::Vector& a0,
                           const Dim<2>::Vector& a1,
                           const Dim<2>::FacetedVolume& poly,
                           const double tol = 1.0e-8);

// 3-D
bool segmentIntersectEdges(const Dim<3>::Vector& a0,
                           const Dim<3>::Vector& a1,
                           const Dim<3>::FacetedVolume& poly,
                           const double tol = 1.0e-8);

}

#endif
