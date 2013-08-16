#include "Utilities/DBC.hh"
#include "Geometry/Dimension.hh"

namespace Spheral {
namespace BoundarySpace {

//------------------------------------------------------------------------------
// Do we have an entry for the NodeList?
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
Boundary<Dimension>::
haveNodeList(const NodeSpace::NodeList<Dimension>& nodeList) const {
  return mBoundaryNodes.find(const_cast<NodeSpace::NodeList<Dimension>*>(&nodeList)) != mBoundaryNodes.end();
}

//------------------------------------------------------------------------------
// Apply the Boundary condtion to the given FieldList.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename DataType>
inline
void
Boundary<Dimension>::
applyFieldListGhostBoundary(FieldSpace::FieldList<Dimension, DataType>& fieldList) const {
  for (typename FieldSpace::FieldList<Dimension, DataType>::iterator fieldItr = fieldList.begin();
       fieldItr < fieldList.end();
       ++fieldItr) {
    REQUIRE(mBoundaryNodes.find(const_cast<NodeSpace::NodeList<Dimension>*>((*fieldItr)->nodeListPtr())) != mBoundaryNodes.end());
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
enforceFieldListBoundary(FieldSpace::FieldList<Dimension, DataType>& fieldList) const {
  for (typename FieldSpace::FieldList<Dimension, DataType>::iterator fieldItr = fieldList.begin();
       fieldItr < fieldList.end();
       ++fieldItr) {
    REQUIRE(mBoundaryNodes.find(const_cast<NodeSpace::NodeList<Dimension>*>((*fieldItr)->nodeListPtr())) != mBoundaryNodes.end());
    enforceBoundary(**fieldItr);
  }
}

//------------------------------------------------------------------------------
// Default choice for whether to mesh ghost nodes.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
Boundary<Dimension>::
meshGhostNodes() const {
  return true;
}

}
}
