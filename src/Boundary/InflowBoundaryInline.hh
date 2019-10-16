#include <typeinfo>

namespace Spheral {

//------------------------------------------------------------------------------
// Return the number of nodes this Boundary is going to create.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
InflowBoundary<Dimension>::
numInflowNodes() const {
  return mNumInflowNodes;
}

//------------------------------------------------------------------------------
// Return the NodeList this Boundary is defined on.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const NodeList<Dimension>&
InflowBoundary<Dimension>::
nodeList() const {
  CHECK(mNodeListPtr != 0);
  return *mNodeListPtr;
}

//------------------------------------------------------------------------------
// The plane for points to enter by.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const GeomPlane<Dimension>& 
InflowBoundary<Dimension>::
plane() const {
  return mPlane;
}

}
