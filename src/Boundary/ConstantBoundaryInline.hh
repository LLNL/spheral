#include <typeinfo>

namespace Spheral {
namespace BoundarySpace {

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

}
}
