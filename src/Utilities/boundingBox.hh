//---------------------------------Spheral++----------------------------------//
// boundingBox
//
// Compute the minimum bounding box for a set of points.
//
// Created by JMO, Wed Oct 19 09:58:58 PDT 2011
//----------------------------------------------------------------------------//
#ifndef __Spheral_boundingBox__
#define __Spheral_boundingBox__

#include <vector>

namespace Spheral {

// Minimum bounding box for a vector of positions.
template<typename Vector>
void
boundingBox(const std::vector<Vector>& positions,
            Vector& xmin,
            Vector& xmax,
            const bool quantize = true);

}

#endif
