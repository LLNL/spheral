//---------------------------------Spheral++----------------------------------//
// ConstantVelocityBoundary -- A boundary condition to enforce a constant 
// velocity on a given set of nodes.
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
#include "Field/FieldList.hh"
#include "Field/Field.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Utilities/DBC.hh"

#include "ConstantVelocityBoundary.hh"

namespace Spheral {
namespace BoundarySpace {

using namespace std;

using NodeSpace::NodeList;
using FieldSpace::Field;
using FieldSpace::FieldList;
using FileIOSpace::FileIO;

//------------------------------------------------------------------------------
// Construct with the given set of nodes.
//------------------------------------------------------------------------------
template<typename Dimension>
ConstantVelocityBoundary<Dimension>::
ConstantVelocityBoundary(const NodeList<Dimension>& nodeList,
                         const vector<int>& nodeIndices):
  Boundary<Dimension>(),
  mNodeListPtr(&nodeList),
  mNodes("Constant Nodes", nodeList, 0),
  mVelocity("Constant velocities", nodeList.velocity()),
  mRestart(DataOutput::registerWithRestart(*this)) {

  // Store the ids of the nodes we're watching.
  for (vector<int>::const_iterator itr = nodeIndices.begin();
       itr < nodeIndices.end();
       ++itr) {
    REQUIRE(*itr >= 0.0 && *itr < nodeList.numInternalNodes());
    mNodes(*itr) = 1;
  }

  // Once we're done the boundary condition should be in a valid state.
  ENSURE(valid());
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
ConstantVelocityBoundary<Dimension>::~ConstantVelocityBoundary() {
}

//------------------------------------------------------------------------------
// Provide the setGhostNodes for a NodeList.  A no-op for this boundary 
// condition.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConstantVelocityBoundary<Dimension>::
setGhostNodes(NodeList<Dimension>& nodeList) {
  this->addNodeList(nodeList);
}

//------------------------------------------------------------------------------
// Provide the updateGhostNodes for a NodeList.  A no-op for this boundary 
// condition.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConstantVelocityBoundary<Dimension>::
updateGhostNodes(NodeList<Dimension>& nodeList) {
}

//------------------------------------------------------------------------------
// Apply the boundary condition to fields of different DataTypes.
// Since this boundary has no ghost nodes, these are all no-ops.
//------------------------------------------------------------------------------
// Specialization for int fields.
template<typename Dimension>
void
ConstantVelocityBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, int>& field) const {
}

// Specialization for scalar fields.
template<typename Dimension>
void
ConstantVelocityBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::Scalar>& field) const {
}

// Specialization for Vector fields.
template<typename Dimension>
void
ConstantVelocityBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::Vector>& field) const {
}

// Specialization for Tensor fields.
template<typename Dimension>
void
ConstantVelocityBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::Tensor>& field) const {
}

// Specialization for symmetric tensors.
template<typename Dimension>
void
ConstantVelocityBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::SymTensor>& field) const {
}

// Specialization for third rank tensors.
template<typename Dimension>
void
ConstantVelocityBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::ThirdRankTensor>& field) const {
}

//------------------------------------------------------------------------------
// Provide the setViolationNodes for a NodeList.  The violation nodes are always
// forced to be the set of nodes we're controlling.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConstantVelocityBoundary<Dimension>::
setViolationNodes(NodeList<Dimension>& nodeList) {

  typedef typename Boundary<Dimension>::BoundaryNodes BoundaryNodes;
  this->addNodeList(nodeList);

  if (&nodeList == mNodeListPtr) {
    BoundaryNodes& boundaryNodes = this->accessBoundaryNodes(nodeList);
    vector<int>& vNodes = boundaryNodes.violationNodes;
    vNodes = nodeIndices();
  }
}

//------------------------------------------------------------------------------
// Provide the updateViolationNodes for a NodeList.  A no-op for this boundary 
// condition.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConstantVelocityBoundary<Dimension>::
updateViolationNodes(NodeList<Dimension>& nodeList) {
}

//------------------------------------------------------------------------------
// Apply the boundary condition to fields of different DataTypes.
//------------------------------------------------------------------------------
// Specialization for int fields, no-op.
template<typename Dimension>
void
ConstantVelocityBoundary<Dimension>::
enforceBoundary(Field<Dimension, int>& field) const {
}

// Specialization for scalar fields, no-op.
template<typename Dimension>
void
ConstantVelocityBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::Scalar>& field) const {
}

// Specialization for Vector fields.  In this case check if we're applying to the
// velocity field or not.  If so, set the velocities to their constant values, otherwise
// a no-op.
template<typename Dimension>
void
ConstantVelocityBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::Vector>& field) const {

  REQUIRE(valid());

  // Is this field the velocity on the NodeList we're watching?
  const NodeList<Dimension>* nodeListPtr = field.nodeListPtr();
  if (nodeListPtr == mNodeListPtr &&
      field.name() == HydroFieldNames::velocity) {

    // This is the velocity field, so enforce the boundary.
    const vector<int> nodeIDs = nodeIndices();
    for (vector<int>::const_iterator itr = nodeIDs.begin();
         itr < nodeIDs.end();
         ++itr) {
      CHECK(*itr < field.numElements());
      field[*itr] = mVelocity[*itr];
    }
  }
}

// Specialization for Tensor fields, no-op.
template<typename Dimension>
void
ConstantVelocityBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::Tensor>& field) const {
}

// Specialization for symmetric tensors, no-op.
template<typename Dimension>
void
ConstantVelocityBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::SymTensor>& field) const {
}

// Specialization for third rank tensors, no-op.
template<typename Dimension>
void
ConstantVelocityBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::ThirdRankTensor>& field) const {
}

//------------------------------------------------------------------------------
// Test if the constant boundary is minimally valid.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
ConstantVelocityBoundary<Dimension>::valid() const {
  return nodeIndices().size() == velocityCondition().size();
}

//------------------------------------------------------------------------------
// Dump the state.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConstantVelocityBoundary<Dimension>::
dumpState(FileIO& file, const string& pathName) const {

  file.write(mNodes, pathName + "/nodeIDs");
  file.write(mVelocity, pathName + "/velocities");

}

//------------------------------------------------------------------------------
// Read the state.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConstantVelocityBoundary<Dimension>::
restoreState(const FileIO& file, const string& pathName) {

  file.read(mNodes, pathName + "/nodeIDs");
  file.read(mVelocity, pathName + "/velocities");

}

}
}
