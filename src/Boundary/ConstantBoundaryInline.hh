#include <typeinfo>

namespace Spheral {
namespace BoundarySpace {

//------------------------------------------------------------------------------
// Return the set of node IDs we're controlling.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const std::vector<int>&
ConstantBoundary<Dimension>::
nodeIndices() const {
  return mNodeIDs;
}

//------------------------------------------------------------------------------
// Return the number of nodes this Boundary is going to create.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
ConstantBoundary<Dimension>::
numConstantNodes() const {
  return mNumConstantNodes;
}

//------------------------------------------------------------------------------
// Return the NodeList this Boundary is defined on.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const NodeSpace::NodeList<Dimension>&
ConstantBoundary<Dimension>::
nodeList() const {
  CHECK(mNodeListPtr != 0);
  return *mNodeListPtr;
}

//------------------------------------------------------------------------------
// The reflection operator.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Tensor
ConstantBoundary<Dimension>::
reflectOperator() const {
  return mReflectOperator;
}

}
}
