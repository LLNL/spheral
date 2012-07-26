//---------------------------------Spheral++----------------------------------//
// ConstantBoundary -- A boundary condition to enforce a constant 
// velocity on a given set of nodes.
//----------------------------------------------------------------------------//
#include "ConstantBoundary.hh"
#include "Boundary.hh"
#include "Field/FieldList.hh"
#include "Field/Field.hh"
#include "Field/FieldBase.hh"

#include "DBC.hh"

namespace Spheral {
namespace BoundarySpace {

using namespace std;
using std::vector;
using std::map;

using NodeSpace::NodeList;
using FieldSpace::FieldBase;
using FieldSpace::Field;
using FieldSpace::FieldList;
using DataBaseSpace::DataBase;

//------------------------------------------------------------------------------
// Construct with the given set of nodes.
//------------------------------------------------------------------------------
template<typename Dimension>
ConstantBoundary<Dimension>::
ConstantBoundary(const NodeList<Dimension>& nodeList,
                 const vector<int>& nodeIDs):
  Boundary<Dimension>(),
  mNumConstantNodes(nodeIDs.size()),
  mNodeListPtr(&nodeList),
  mScalarValues(),
  mVectorValues(),
  mVector3dValues(),
  mTensorValues(),
  mSymTensorValues(),
  mThirdRankTensorValues(),
  mVectorScalarValues() {

  // Store the field values on the given nodes.
  storeFieldValues<int>(nodeList, nodeIDs, mIntValues);
  storeFieldValues<Scalar>(nodeList, nodeIDs, mScalarValues);
  storeFieldValues<Vector>(nodeList, nodeIDs, mVectorValues);
  storeFieldValues<Vector3d>(nodeList, nodeIDs, mVector3dValues);
  storeFieldValues<Tensor>(nodeList, nodeIDs, mTensorValues);
  storeFieldValues<SymTensor>(nodeList, nodeIDs, mSymTensorValues);
  storeFieldValues<ThirdRankTensor>(nodeList, nodeIDs, mThirdRankTensorValues);
  storeFieldValues<std::vector<Scalar> >(nodeList, nodeIDs, mVectorScalarValues);

  // Once we're done the boundary condition should be in a valid state.
  ENSURE(valid());
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
ConstantBoundary<Dimension>::~ConstantBoundary() {
}

//------------------------------------------------------------------------------
// Provide the setGhostNodes for a NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConstantBoundary<Dimension>::
setGhostNodes(NodeList<Dimension>& nodeList) {

  // Add the control/ghost nodes for this NodeList.
  this->addNodeList(nodeList);
  typename Boundary<Dimension>::BoundaryNodes& boundaryNodes = this->accessBoundaryNodes(nodeList);
  vector<int>& ghostNodes = boundaryNodes.ghostNodes;
  ghostNodes.resize(0);

  // Check to see if this is the NodeList we're working on.
  if (&nodeList == mNodeListPtr) {

    // Set the ghost node ID's we'll be controlling.
    CHECK(ghostNodes.size() == 0);
    ghostNodes.reserve(mNumConstantNodes);
    const int currentNumGhostNodes = nodeList.numGhostNodes();
    const int firstNewGhostNode = nodeList.numNodes();
    nodeList.numGhostNodes(currentNumGhostNodes + mNumConstantNodes);
    CHECK(nodeList.numNodes() == firstNewGhostNode + mNumConstantNodes);
    for (int i = 0; i < mNumConstantNodes; ++i) {
      ghostNodes.push_back(firstNewGhostNode + i);
    }
    CHECK(ghostNodes.size() == mNumConstantNodes);

    // Set us up as our own control nodes.
    boundaryNodes.controlNodes = ghostNodes;

    // Set the ghost node positions and H.
    this->updateGhostNodes(nodeList);
  }
}

//------------------------------------------------------------------------------
// Provide the updateGhostNodes for a NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConstantBoundary<Dimension>::
updateGhostNodes(NodeList<Dimension>& nodeList) {

  // Check to see if this is the NodeList we're working on.
  if (&nodeList == mNodeListPtr) {

    // Set the ghost node positions.
    Field<Dimension, Vector>& positions = nodeList.positions();
    setGhostValues<Vector>(positions, mVectorValues);

    // Set the ghost H fields.
    Field<Dimension, SymTensor>& Hfield = nodeList.Hfield();
    setGhostValues<SymTensor>(Hfield, mSymTensorValues);

//     // Update the neighbor information.
//     const vector<int>& ghostNodes = this->ghostNodes(nodeList);
//     CHECK(ghostNodes.size() == mNumConstantNodes);
//     nodeList.neighbor().updateNodes(); // (ghostNodes);
  }
}

//------------------------------------------------------------------------------
// Apply the boundary condition to fields of different DataTypes.
//------------------------------------------------------------------------------
// Specialization for int fields.
template<typename Dimension>
void
ConstantBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, int>& field) const {
  setGhostValues<int>(field, mIntValues);
}

// Specialization for scalar fields.
template<typename Dimension>
void
ConstantBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::Scalar>& field) const {
  setGhostValues<Scalar>(field, mScalarValues);
}

// Specialization for Vector fields.
template<typename Dimension>
void
ConstantBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::Vector>& field) const {
  setGhostValues<Vector>(field, mVectorValues);
}

// Specialization for Vector3d fields.
template<typename Dimension>
void
ConstantBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::Vector3d>& field) const {
  setGhostValues<Vector3d>(field, mVector3dValues);
}

// Specialization for Tensor fields.
template<typename Dimension>
void
ConstantBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::Tensor>& field) const {
  setGhostValues<Tensor>(field, mTensorValues);
}

// Specialization for symmetric tensors.
template<typename Dimension>
void
ConstantBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::SymTensor>& field) const {
  setGhostValues<SymTensor>(field, mSymTensorValues);
}

// Specialization for third rank tensors.
template<typename Dimension>
void
ConstantBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::ThirdRankTensor>& field) const {
  setGhostValues<ThirdRankTensor>(field, mThirdRankTensorValues);
}

// Specialization for vector<scalar> fields.
template<typename Dimension>
void
ConstantBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, std::vector<typename Dimension::Scalar> >& field) const {
  setGhostValues<std::vector<Scalar> >(field, mVectorScalarValues);
}

//------------------------------------------------------------------------------
// Provide the setViolationNodes for a NodeList.  A no-op for this boundary 
// condition.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConstantBoundary<Dimension>::
setViolationNodes(NodeList<Dimension>& nodeList) {
}

//------------------------------------------------------------------------------
// Provide the updateViolationNodes for a NodeList.  A no-op for this boundary 
// condition.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConstantBoundary<Dimension>::
updateViolationNodes(NodeList<Dimension>& nodeList) {
}

//------------------------------------------------------------------------------
// Apply the boundary condition to fields of different DataTypes.  A no-op.
//------------------------------------------------------------------------------
// Specialization for int fields.
template<typename Dimension>
void
ConstantBoundary<Dimension>::
enforceBoundary(Field<Dimension, int>& field) const {
}

// Specialization for scalar fields.
template<typename Dimension>
void
ConstantBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::Scalar>& field) const {
}

// Specialization for Vector fields.
template<typename Dimension>
void
ConstantBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::Vector>& field) const {
}

// Specialization for Vector fields.
template<typename Dimension>
void
ConstantBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::Vector3d>& field) const {
}

// Specialization for Tensor fields.
template<typename Dimension>
void
ConstantBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::Tensor>& field) const {
}

// Specialization for symmetric tensors.
template<typename Dimension>
void
ConstantBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::SymTensor>& field) const {
}

// Specialization for third rank tensors.
template<typename Dimension>
void
ConstantBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::ThirdRankTensor>& field) const {
}

//------------------------------------------------------------------------------
// Test if the constant boundary is minimally valid.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
ConstantBoundary<Dimension>::valid() const {
  bool ok = true;
  for (typename std::map<const FieldBase<Dimension>*, std::vector<int> >::const_iterator
         itr = mIntValues.begin();
       itr != mIntValues.end();
       ++itr) ok = ok && itr->first->nodeListPtr() == mNodeListPtr;
  for (typename std::map<const FieldBase<Dimension>*, std::vector<Scalar> >::const_iterator
         itr = mScalarValues.begin();
       itr != mScalarValues.end();
       ++itr) ok = ok && itr->first->nodeListPtr() == mNodeListPtr;
  for (typename std::map<const FieldBase<Dimension>*, std::vector<Vector> >::const_iterator
         itr = mVectorValues.begin();
       itr != mVectorValues.end();
       ++itr) ok = ok && itr->first->nodeListPtr() == mNodeListPtr;
  for (typename std::map<const FieldBase<Dimension>*, std::vector<Tensor> >::const_iterator
         itr = mTensorValues.begin();
       itr != mTensorValues.end();
       ++itr) ok = ok && itr->first->nodeListPtr() == mNodeListPtr;
  for (typename std::map<const FieldBase<Dimension>*, std::vector<SymTensor> >::const_iterator
         itr = mSymTensorValues.begin();
       itr != mSymTensorValues.end();
       ++itr) ok = ok && itr->first->nodeListPtr() == mNodeListPtr;
  for (typename std::map<const FieldBase<Dimension>*, std::vector<ThirdRankTensor> >::const_iterator
         itr = mThirdRankTensorValues.begin();
       itr != mThirdRankTensorValues.end();
       ++itr) ok = ok && itr->first->nodeListPtr() == mNodeListPtr;
  for (typename std::map<const FieldBase<Dimension>*, std::vector<std::vector<Scalar> > >::const_iterator
         itr = mVectorScalarValues.begin();
       itr != mVectorScalarValues.end();
       ++itr) ok = ok && itr->first->nodeListPtr() == mNodeListPtr;
  return ok;
}

}
}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"

namespace Spheral {
namespace BoundarySpace {
template class ConstantBoundary< Dim<1> >;
template class ConstantBoundary< Dim<2> >;
template class ConstantBoundary< Dim<3> >;
}
}
