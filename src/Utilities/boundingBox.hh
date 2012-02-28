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

// Forward declarations.
namespace Spheral {
  namespace FieldSpace {
    template<typename Dimension, typename Value> class FieldList;
  }
}

namespace Spheral {

// Minimum bounding box for a vector of positions.
template<typename Vector>
void
boundingBox(const std::vector<Vector>& positions,
            Vector& xmin,
            Vector& xmax,
            const bool quantize = true);

template<typename Dimension>
void
boundingBox(const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& positions,
            typename Dimension::Vector& xmin,
            typename Dimension::Vector& xmax,
            const bool quantize,
            const bool useGhosts);

}

#endif
