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

template<typename Dimension, typename Value> class FieldList;

// Minimum bounding box for a vector of positions.
template<typename Vector>
void
boundingBox(const std::vector<Vector>& positions,
            Vector& xmin,
            Vector& xmax);

template<typename Dimension>
void
boundingBox(const FieldList<Dimension, typename Dimension::Vector>& positions,
            typename Dimension::Vector& xmin,
            typename Dimension::Vector& xmax,
            const bool useGhosts);

}

#endif
