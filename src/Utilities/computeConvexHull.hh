//---------------------------------Spheral++----------------------------------//
// computeConvexHull
//
// Returns the convex hull for a set of points.
//
// Created by JMO, Tue Jan 26 14:33:35 PST 2010
//----------------------------------------------------------------------------//
#ifndef __Spheral_computeConvexHull__
#define __Spheral_computeConvexHull__

#include <vector>
#include "Geometry/Dimension.hh"

namespace Spheral {

Dim<1>::ConvexHull
computeConvexHull(const std::vector<Dim<1>::WMVector>& points);

Dim<2>::ConvexHull
computeConvexHull(const std::vector<Dim<2>::WMVector>& points);

Dim<3>::ConvexHull
computeConvexHull(const std::vector<Dim<3>::WMVector>& points);

}

#endif
