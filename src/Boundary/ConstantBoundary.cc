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
#include "cdebug.hh"

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
  cdebug << "ConstantBoundary::ConstantBoundary" << this << endl;

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
  cdebug << "ConstantBoundary::~ConstantBoundary() " << this << endl;
}

//------------------------------------------------------------------------------
// Provide the setGhostNodes for a NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConstantBoundary<Dimension>::
setGhostNodes(NodeList<Dimension>& nodeList) {
  cdebug << "ConstantBoundary::setGhostNodes(NodeList<Dimension>&)" << endl;

  // Add the control/ghost nodes for this NodeList.
  addNodeList(nodeList);
  typename Boundary<Dimension>::BoundaryNodes& boundaryNodes = accessBoundaryNodes(nodeList);
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

    // Set the ghost node positions and H.
    updateGhostNodes(nodeList);
  }
}

//------------------------------------------------------------------------------
// Provide the updateGhostNodes for a NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConstantBoundary<Dimension>::
updateGhostNodes(NodeList<Dimension>& nodeList) {
  cdebug << "ConstantBoundary::updateGhostNodes(NodeList<Dimension>&)" << endl;

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
  cdebug << "ConstantBoundary::applyGhostBoundary(IntField) " << this << endl;
  setGhostValues<int>(field, mIntValues);
}

// Specialization for scalar fields.
template<typename Dimension>
void
ConstantBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::Scalar>& field) const {
  cdebug << "ConstantBoundary::applyGhostBoundary(ScalarField) " << this << endl;
  setGhostValues<Scalar>(field, mScalarValues);
}

// Specialization for Vector fields.
template<typename Dimension>
void
ConstantBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::Vector>& field) const {
  cdebug << "ConstantBoundary::applyGhostBoundary(VectorField) " << this << endl;
  setGhostValues<Vector>(field, mVectorValues);
}

// Specialization for Vector3d fields.
template<typename Dimension>
void
ConstantBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::Vector3d>& field) const {
  cdebug << "ConstantBoundary::applyVector3dGhostBoundary(VectorField) " << this << endl;
  setGhostValues<Vector3d>(field, mVector3dValues);
}

// Specialization for Tensor fields.
template<typename Dimension>
void
ConstantBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::Tensor>& field) const {
  cdebug << "ConstantBoundary::applyGhostBoundary(VectorField) " << this << endl;
  setGhostValues<Tensor>(field, mTensorValues);
}

// Specialization for symmetric tensors.
template<typename Dimension>
void
ConstantBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::SymTensor>& field) const {
  cdebug << "ConstantBoundary::applyGhostBoundary(SymTensorField) " << this << endl;
  setGhostValues<SymTensor>(field, mSymTensorValues);
}

// Specialization for third rank tensors.
template<typename Dimension>
void
ConstantBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::ThirdRankTensor>& field) const {
  cdebug << "ConstantBoundary::applyGhostBoundary(ThirdRankTensorField) " << this << endl;
  setGhostValues<ThirdRankTensor>(field, mThirdRankTensorValues);
}

// Specialization for vector<scalar> fields.
template<typename Dimension>
void
ConstantBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, std::vector<typename Dimension::Scalar> >& field) const {
  cdebug << "ConstantBoundary::applyGhostBoundary(VectorScalarField) " << this << endl;
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
  cdebug << "ConstantBoundary::setViolationNodes(NodeList<Dimension>&)" << endl;
}

//------------------------------------------------------------------------------
// Provide the updateViolationNodes for a NodeList.  A no-op for this boundary 
// condition.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConstantBoundary<Dimension>::
updateViolationNodes(NodeList<Dimension>& nodeList) {
  cdebug << "ConstantBoundary::updateViolationNodes(NodeList<Dimension>&)" << endl;
}

//------------------------------------------------------------------------------
// Apply the boundary condition to fields of different DataTypes.  A no-op.
//------------------------------------------------------------------------------
// Specialization for int fields.
template<typename Dimension>
void
ConstantBoundary<Dimension>::
enforceBoundary(Field<Dimension, int>& field) const {
  cdebug << "ConstantBoundary::enforceBoundary(IntField) " << this << endl;
}

// Specialization for scalar fields.
template<typename Dimension>
void
ConstantBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::Scalar>& field) const {
  cdebug << "ConstantBoundary::enforceBoundary(ScalarField) " << this << endl;
}

// Specialization for Vector fields.
template<typename Dimension>
void
ConstantBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::Vector>& field) const {
  cdebug << "ConstantBoundary::enforceBoundary(VectorField) " << this << endl;
}

// Specialization for Vector fields.
template<typename Dimension>
void
ConstantBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::Vector3d>& field) const {
  cdebug << "ConstantBoundary::enforceVector3dBoundary(Vector3dField) " << this << endl;
}

// Specialization for Tensor fields.
template<typename Dimension>
void
ConstantBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::Tensor>& field) const {
  cdebug << "ConstantBoundary::enforceBoundary(VectorField) " << this << endl;
}

// Specialization for symmetric tensors.
template<typename Dimension>
void
ConstantBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::SymTensor>& field) const {
  cdebug << "ConstantBoundary::enforceBoundary(SymTensorField) " << this << endl;
}

// Specialization for third rank tensors.
template<typename Dimension>
void
ConstantBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::ThirdRankTensor>& field) const {
  cdebug << "ConstantBoundary::enforceBoundary(ThirdRankTensorField) " << this << endl;
}

//------------------------------------------------------------------------------
// Test if the constant boundary is minimally valid.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
ConstantBoundary<Dimension>::valid() const {
  cdebug << "ConstantBoundary::valid() " << this << endl;
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
