#include "Geometry/Dimension.hh"
#include <vector>

namespace Spheral {

//------------------------------------------------------------------------------
// Find the the matching vertex between two lists (assumed to be in opposite 
// order.)
//------------------------------------------------------------------------------
unsigned
findMatchingVertex(const unsigned i, 
                   const std::vector<Dim<2>::Vector>& verticesi,
                   const std::vector<Dim<2>::Vector>& verticesj);


//------------------------------------------------------------------------------
// Find the vertex in the second list that best lines up with the first in 
// the first list.  Assumes that the two lists are the same vertices in opposite
// order with arbitrary starting points.
//------------------------------------------------------------------------------
unsigned
findMatchingVertex(const std::vector<Dim<3>::Vector>& verticesi,
                   const std::vector<Dim<3>::Vector>& verticesj,
                   const std::vector<unsigned>& indicesi,
                   const std::vector<unsigned>& indicesj);

//------------------------------------------------------------------------------
// Find the closest vertex in the list to the given position.
//------------------------------------------------------------------------------
template<typename Vector>
unsigned
findMatchingVertex(const Vector& target,
                   const std::vector<Vector>& verticesj);

//------------------------------------------------------------------------------
// Find the closest vertex in a subset of the list to the given position.
//------------------------------------------------------------------------------
template<typename Vector>
unsigned
findMatchingVertex(const Vector& target,
                   const std::vector<Vector>& verticesj,
                   const std::vector<unsigned>& indicesj);

}
