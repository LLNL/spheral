//---------------------------------Spheral++----------------------------------//
// ConstantVelocityBoundary -- A boundary condition to enforce a constant 
// velocity on a given set of nodes.
//----------------------------------------------------------------------------//
#include "ConstantVelocityBoundary.hh"
#include "Field/FieldList.hh"
#include "Field/Field.hh"
#include "Hydro/HydroFieldNames.hh"
#include "FileIO/FileIO.hh"

#include "DBC.hh"
#include "cdebug.hh"

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
                         const vector<int>& nodeIndicies):
  Boundary<Dimension>(),
  mNodeListPtr(&nodeList),
  mNodes("Constant Nodes", nodeList, 0),
  mVelocity("Constant velocities", nodeList.velocity()),
  mRestart(DataOutput::registerWithRestart(*this)) {
  cdebug << "ConstantVelocityBoundary::ConstantVelocityBoundary" << this << endl;

  // Store the ids of the nodes we're watching.
  for (vector<int>::const_iterator itr = nodeIndicies.begin();
       itr < nodeIndicies.end();
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
  cdebug << "ConstantVelocityBoundary::~ConstantVelocityBoundary() " << this << endl;
}

//------------------------------------------------------------------------------
// Provide the setGhostNodes for a NodeList.  A no-op for this boundary 
// condition.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConstantVelocityBoundary<Dimension>::
setGhostNodes(NodeList<Dimension>& nodeList) {
  cdebug << "ConstantVelocityBoundary::setGhostNodes(NodeList<Dimension>&)" << endl;
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
  cdebug << "ConstantVelocityBoundary::updateGhostNodes(NodeList<Dimension>&)" << endl;
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
  cdebug << "ConstantVelocityBoundary::applyGhostBoundary(IntField) " << this << endl;
}

// Specialization for scalar fields.
template<typename Dimension>
void
ConstantVelocityBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::Scalar>& field) const {
  cdebug << "ConstantVelocityBoundary::applyGhostBoundary(ScalarField) " << this << endl;
}

// Specialization for Vector fields.
template<typename Dimension>
void
ConstantVelocityBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::Vector>& field) const {
  cdebug << "ConstantVelocityBoundary::applyGhostBoundary(VectorField) " << this << endl;
}

// Specialization for Vector3d fields.
template<typename Dimension>
void
ConstantVelocityBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::Vector3d>& field) const {
  cdebug << "ConstantVelocityBoundary::applyGhostBoundary(Vector3dField) " << this << endl;
}

// Specialization for Tensor fields.
template<typename Dimension>
void
ConstantVelocityBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::Tensor>& field) const {
  cdebug << "ConstantVelocityBoundary::applyGhostBoundary(VectorField) " << this << endl;
}

// Specialization for symmetric tensors.
template<typename Dimension>
void
ConstantVelocityBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::SymTensor>& field) const {
  cdebug << "ConstantVelocityBoundary::applyGhostBoundary(SymTensorField) " << this << endl;
}

// Specialization for third rank tensors.
template<typename Dimension>
void
ConstantVelocityBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::ThirdRankTensor>& field) const {
  cdebug << "ConstantVelocityBoundary::applyGhostBoundary(ThirdRankTensorField) " << this << endl;
}

// Specialization for vector<scalar> fields.
template<typename Dimension>
void
ConstantVelocityBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, std::vector<typename Dimension::Scalar> >& field) const {
  cdebug << "ConstantVelocityBoundary::applyGhostBoundary(VectorScalarField) " << this << endl;
}

//------------------------------------------------------------------------------
// Provide the setViolationNodes for a NodeList.  The violation nodes are always
// forced to be the set of nodes we're controlling.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConstantVelocityBoundary<Dimension>::
setViolationNodes(NodeList<Dimension>& nodeList) {
  cdebug << "ConstantVelocityBoundary::setViolationNodes(NodeList<Dimension>&)" << endl;

  typedef typename Boundary<Dimension>::BoundaryNodes BoundaryNodes;
  this->addNodeList(nodeList);

  if (&nodeList == mNodeListPtr) {
    BoundaryNodes& boundaryNodes = this->accessBoundaryNodes(nodeList);
    vector<int>& vNodes = boundaryNodes.violationNodes;
    vNodes = nodeIndicies();
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
  cdebug << "ConstantVelocityBoundary::updateViolationNodes(NodeList<Dimension>&)" << endl;
}

//------------------------------------------------------------------------------
// Apply the boundary condition to fields of different DataTypes.
//------------------------------------------------------------------------------
// Specialization for int fields, no-op.
template<typename Dimension>
void
ConstantVelocityBoundary<Dimension>::
enforceBoundary(Field<Dimension, int>& field) const {
  cdebug << "ConstantVelocityBoundary::enforceBoundary(IntField) " << this << endl;
}

// Specialization for scalar fields, no-op.
template<typename Dimension>
void
ConstantVelocityBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::Scalar>& field) const {
  cdebug << "ConstantVelocityBoundary::enforceBoundary(ScalarField) " << this << endl;
}

// Specialization for Vector fields.  In this case check if we're applying to the
// velocity field or not.  If so, set the velocities to their constant values, otherwise
// a no-op.
template<typename Dimension>
void
ConstantVelocityBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::Vector>& field) const {
  cdebug << "ConstantVelocityBoundary::enforceBoundary(VectorField) " << this << endl;

  REQUIRE(valid());

  // Is this field the velocity on the NodeList we're watching?
  const NodeList<Dimension>* nodeListPtr = field.nodeListPtr();
  if (nodeListPtr == mNodeListPtr &&
      field.name() == HydroFieldNames::velocity) {

    // This is the velocity field, so enforce the boundary.
    const vector<int> nodeIDs = nodeIndicies();
    for (vector<int>::const_iterator itr = nodeIDs.begin();
         itr < nodeIDs.end();
         ++itr) {
      CHECK(*itr < field.numElements());
      field[*itr] = mVelocity[*itr];
    }
  }
}

// Specialization for vector3d fields, no-op.
template<typename Dimension>
void
ConstantVelocityBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::Vector3d>& field) const {
  cdebug << "ConstantVelocityBoundary::enforceBoundary(Vector3dField) " << this << endl;
}

// Specialization for Tensor fields, no-op.
template<typename Dimension>
void
ConstantVelocityBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::Tensor>& field) const {
  cdebug << "ConstantVelocityBoundary::enforceBoundary(TensorField) " << this << endl;
}

// Specialization for symmetric tensors, no-op.
template<typename Dimension>
void
ConstantVelocityBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::SymTensor>& field) const {
  cdebug << "ConstantVelocityBoundary::enforceBoundary(SymTensorField) " << this << endl;
}

// Specialization for third rank tensors, no-op.
template<typename Dimension>
void
ConstantVelocityBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::ThirdRankTensor>& field) const {
  cdebug << "ConstantVelocityBoundary::enforceBoundary(ThirdRankTensorField) " << this << endl;
}

//------------------------------------------------------------------------------
// Test if the constant boundary is minimally valid.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
ConstantVelocityBoundary<Dimension>::valid() const {
  cdebug << "ConstantVelocityBoundary::valid() " << this << endl;
  return nodeIndicies().size() == velocityCondition().size();
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
