#include <vector>
#include "Geometry/Dimension.hh"

namespace Spheral {
namespace MeshSpace {

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
unsigned
findMatchingVertex(const Dim<3>::Vector& target,
                   const std::vector<Dim<3>::Vector>& verticesj);

//------------------------------------------------------------------------------
// Find the closest vertex in a subset of the list to the given position.
//------------------------------------------------------------------------------
unsigned
findMatchingVertex(const Dim<3>::Vector& target,
                   const std::vector<Dim<3>::Vector>& verticesj,
                   const std::vector<unsigned>& indicesj);

}
}
