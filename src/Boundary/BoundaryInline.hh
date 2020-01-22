#include "Geometry/Dimension.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Apply the Boundary condtion to the given FieldList.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename DataType>
inline
void
Boundary<Dimension>::
applyFieldListGhostBoundary(FieldList<Dimension, DataType>& fieldList) const {
  for (typename FieldList<Dimension, DataType>::iterator fieldItr = fieldList.begin();
       fieldItr < fieldList.end();
       ++fieldItr) {
    REQUIRE(mBoundaryNodes.find(const_cast<NodeList<Dimension>*>((*fieldItr)->nodeListPtr())) != mBoundaryNodes.end());
    applyGhostBoundary(**fieldItr);
  }
}

//------------------------------------------------------------------------------
// Enforce the Boundary condtion on the violation nodes in the given FieldList.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename DataType>
inline
void
Boundary<Dimension>::
enforceFieldListBoundary(FieldList<Dimension, DataType>& fieldList) const {
  for (typename FieldList<Dimension, DataType>::iterator fieldItr = fieldList.begin();
       fieldItr < fieldList.end();
       ++fieldItr) {
    REQUIRE(mBoundaryNodes.find(const_cast<NodeList<Dimension>*>((*fieldItr)->nodeListPtr())) != mBoundaryNodes.end());
    enforceBoundary(**fieldItr);
  }
}

//------------------------------------------------------------------------------
// Do we have an entry for the NodeList?
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
Boundary<Dimension>::
haveNodeList(const NodeList<Dimension>& nodeList) const {
  return mBoundaryNodes.find(const_cast<NodeList<Dimension>*>(&nodeList)) != mBoundaryNodes.end();
}

}
