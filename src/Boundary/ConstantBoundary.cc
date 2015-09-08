//---------------------------------Spheral++----------------------------------//
// ConstantBoundary -- A boundary condition to enforce a constant 
// velocity on a given set of nodes.
//----------------------------------------------------------------------------//
#include "ConstantBoundary.hh"
#include "Boundary.hh"
#include "mapPositionThroughPlanes.hh"
#include "Field/FieldList.hh"
#include "Field/Field.hh"
#include "Field/FieldBase.hh"
#include "Utilities/planarReflectingOperator.hh"
#include "Utilities/DBC.hh"

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

namespace {

//------------------------------------------------------------------------------
// Store the Field values of the given DataType for the given NodeList.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
void
storeFieldValues(const NodeSpace::NodeList<Dimension>& nodeList,
                 const std::vector<int>& nodeIDs,
                 std::map<std::string, std::vector<DataType> >& values) {

  // Iterate over all the Fields defined on the NodeList.
  for (typename NodeSpace::NodeList<Dimension>::const_FieldBaseIterator fieldItr = nodeList.registeredFieldsBegin();
       fieldItr != nodeList.registeredFieldsEnd();
       ++fieldItr) {

    // Determine if this Field is the type we're looking for.
    if (typeid(**fieldItr) == typeid(FieldSpace::Field<Dimension, DataType>)) {
      const FieldSpace::Field<Dimension, DataType>& field = (const FieldSpace::Field<Dimension, DataType>&) **fieldItr;

      // Build a vector of the values of this field on the requested nodes.
      std::vector<DataType> vals;
      vals.reserve(nodeIDs.size());
      for (typename std::vector<int>::const_iterator nodeItr = nodeIDs.begin();
           nodeItr != nodeIDs.end();
           ++nodeItr) {
        CHECK(*nodeItr >= 0 && *nodeItr < field.numElements());
        vals.push_back(field(*nodeItr));
      }
      CHECK(vals.size() == nodeIDs.size());

      // Put this data into the result.
      // cerr << " --> ";
      // for (typename std::map<std::string, std::vector<DataType> >::const_iterator itr = values.begin();
      //      itr != values.end();
      //      ++itr) cerr << itr->first << " ";
      // cerr << endl;
      const std::string key = StateBase<Dimension>::key(**fieldItr);
      CHECK2(values.find(key) == values.end(), key);
      values[key] = vals;
    }
  }
}
  
//------------------------------------------------------------------------------
// Set the Field values in the given field using the given map.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
void
resetValues(FieldSpace::Field<Dimension, DataType>& field,
            const std::vector<int>& nodeIDs,
            const std::map<std::string, std::vector<DataType> >& values,
            const bool dieOnMissingField) {

  const NodeSpace::NodeList<Dimension>& nodeList = field.nodeList();

  // Find this Field in the set of stored values.
  const std::string key = StateBase<Dimension>::key(field);
  typename std::map<std::string, std::vector<DataType> >::const_iterator itr = values.find(key);
  VERIFY2(itr != values.end() or not dieOnMissingField,
          "ConstantBoundary error: " << key << " not found in stored field values.");

  // Now set the values.
  if (itr != values.end()) {
    const std::vector<DataType>& newValues = itr->second;;
    CHECK2(newValues.size() == nodeIDs.size(), newValues.size() << " " << nodeIDs.size());
    for (int i = 0; i < nodeIDs.size(); ++i) {
      CHECK(nodeIDs[i] >= 0 &&
            nodeIDs[i] < nodeList.numNodes());
      field(nodeIDs[i]) = newValues[i];
    }
  }
}

}

//------------------------------------------------------------------------------
// Construct with the given set of nodes.
//------------------------------------------------------------------------------
template<typename Dimension>
ConstantBoundary<Dimension>::
ConstantBoundary(NodeList<Dimension>& nodeList,
                 const vector<int>& nodeIDs,
                 const GeomPlane<Dimension>& denialPlane):
  Boundary<Dimension>(),
  mNodeListPtr(&nodeList),
  mNodeIDs(nodeIDs),
  mNumConstantNodes(nodeIDs.size()),
  mDenialPlane(denialPlane),
  mReflectOperator(planarReflectingOperator(denialPlane)),
  mActive(false),
  mScalarValues(),
  mVectorValues(),
  mTensorValues(),
  mSymTensorValues(),
  mThirdRankTensorValues(),
  mVectorScalarValues() {

  // // Store the ids of the nodes we're watching.
  // for (vector<int>::const_iterator itr = nodeIDs.begin();
  //      itr < nodeIDs.end();
  //      ++itr) {
  //   REQUIRE(*itr >= 0.0 && *itr < nodeList.numInternalNodes());
  //   mNodeFlags[*itr] = 1;
  // }

  // Issue a big old warning!
  if (Process::getRank() == 0) cerr << "WARNING: ConstantBoundary is currently not compatible with redistributing nodes!\nMake sure you don't allow redistribution with this Boundary condition." << endl;

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
// setGhostNodes
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConstantBoundary<Dimension>::
setGhostNodes(NodeList<Dimension>& nodeList) {
  this->addNodeList(nodeList);
  cerr << "1" << endl;

  if (mActive and &nodeList == mNodeListPtr) {
    typename Boundary<Dimension>::BoundaryNodes& boundaryNodes = this->accessBoundaryNodes(nodeList);
    vector<int>& cNodes = boundaryNodes.controlNodes;
    vector<int>& gNodes = boundaryNodes.ghostNodes;
    cerr << "2" << endl;
    
    unsigned currentNumGhostNodes = nodeList.numGhostNodes();
    unsigned firstNewGhostNode = nodeList.numNodes();
    nodeList.numGhostNodes(currentNumGhostNodes + mNumConstantNodes);
    mNodeIDs = vector<int>(mNumConstantNodes);
    cNodes = vector<int>(mNumConstantNodes);
    gNodes = vector<int>(mNumConstantNodes);
    cerr << "3" << endl;
    for (int i = 0; i != mNumConstantNodes; ++i) {
      mNodeIDs[i] = firstNewGhostNode + i;
      cNodes[i] = mNodeIDs[i];
      gNodes[i] = mNodeIDs[i];
    }
    cerr << "4" << endl;
    
    this->updateGhostNodes(nodeList);
    cerr << "5" << endl;
  }
}

//------------------------------------------------------------------------------
// updateGhostNodes
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConstantBoundary<Dimension>::
updateGhostNodes(NodeList<Dimension>& nodeList) {
  if (mActive and &nodeList == mNodeListPtr) {
    this->applyGhostBoundary(nodeList.positions());
    this->applyGhostBoundary(nodeList.Hfield());
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
  if (mActive) resetValues(field, this->nodeIndices(), mIntValues, false);
}

// Specialization for scalar fields.
template<typename Dimension>
void
ConstantBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::Scalar>& field) const {
  if (mActive) resetValues(field, this->nodeIndices(), mScalarValues, false);
}

// Specialization for Vector fields.
template<typename Dimension>
void
ConstantBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::Vector>& field) const {
  if (mActive) resetValues(field, this->nodeIndices(), mVectorValues, false);
}

// Specialization for Tensor fields.
template<typename Dimension>
void
ConstantBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::Tensor>& field) const {
  if (mActive) resetValues(field, this->nodeIndices(), mTensorValues, false);
}

// Specialization for symmetric tensors.
template<typename Dimension>
void
ConstantBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::SymTensor>& field) const {
  if (mActive) resetValues(field, this->nodeIndices(), mSymTensorValues, false);
}

// Specialization for third rank tensors.
template<typename Dimension>
void
ConstantBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::ThirdRankTensor>& field) const {
  if (mActive) resetValues(field, this->nodeIndices(), mThirdRankTensorValues, false);
}

// Specialization for vector<scalar> fields.
template<typename Dimension>
void
ConstantBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, std::vector<typename Dimension::Scalar> >& field) const {
  if (mActive) resetValues(field, this->nodeIndices(), mVectorScalarValues, false);
}

//------------------------------------------------------------------------------
// Provide the setViolationNodes for a NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConstantBoundary<Dimension>::
setViolationNodes(NodeList<Dimension>& nodeList) {

  typedef typename Boundary<Dimension>::BoundaryNodes BoundaryNodes;
  this->addNodeList(nodeList);

  if (&nodeList == mNodeListPtr) {
    BoundaryNodes& boundaryNodes = this->accessBoundaryNodes(nodeList);
    vector<int>& vNodes = boundaryNodes.violationNodes;
    const Field<Dimension, Vector>& pos = nodeList.positions();
    for (unsigned i = 0; i != nodeList.numInternalNodes(); ++i) {
      if (mDenialPlane.compare(pos(i)) == -1) vNodes.push_back(i);
      // if (mNodeFlags[i] == 1 or mDenialPlane.compare(pos(i)) == -1) vNodes.push_back(i);
    }

    this->updateViolationNodes(nodeList);
  }
}

//------------------------------------------------------------------------------
// Provide the updateViolationNodes for a NodeList, correcting positions and 
// H's.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConstantBoundary<Dimension>::
updateViolationNodes(NodeList<Dimension>& nodeList) {
  if (&nodeList == mNodeListPtr) {
    Field<Dimension, Vector>& pos = nodeList.positions();
    Field<Dimension, Vector>& vel = nodeList.velocity();
    Field<Dimension, SymTensor>& H = nodeList.Hfield();
    // this->enforceBoundary(pos);
    // this->enforceBoundary(H);

    // Look for any nodes in violation of the plane and reset their positions
    // and H's.
    for (unsigned i = 0; i != nodeList.numInternalNodes(); ++i) {
      if (mDenialPlane.compare(pos(i)) == 1) {
        pos(i) = mapPositionThroughPlanes(pos(i), mDenialPlane, mDenialPlane);
        vel(i) = mReflectOperator*vel(i);
        H(i) = (mReflectOperator*H(i)*mReflectOperator).Symmetric();
      }
    }
  }
}

//------------------------------------------------------------------------------
// Apply the boundary condition to fields of different DataTypes.
//------------------------------------------------------------------------------
// Specialization for int fields.
template<typename Dimension>
void
ConstantBoundary<Dimension>::
enforceBoundary(Field<Dimension, int>& field) const {
  // resetValues(field, this->nodeIndices(), mIntValues);
}

// Specialization for scalar fields.
template<typename Dimension>
void
ConstantBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::Scalar>& field) const {
  // resetValues(field, this->nodeIndices(), mScalarValues);
}

// Specialization for Vector fields.
template<typename Dimension>
void
ConstantBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::Vector>& field) const {
  // resetValues(field, this->nodeIndices(), mVectorValues);
}

// Specialization for Tensor fields.
template<typename Dimension>
void
ConstantBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::Tensor>& field) const {
  // resetValues(field, this->nodeIndices(), mTensorValues);
}

// Specialization for symmetric tensors.
template<typename Dimension>
void
ConstantBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::SymTensor>& field) const {
  // resetValues(field, this->nodeIndices(), mSymTensorValues);
}

// Specialization for third rank tensors.
template<typename Dimension>
void
ConstantBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::ThirdRankTensor>& field) const {
  // resetValues(field, this->nodeIndices(), mThirdRankTensorValues);
}

// // Specialization for vector<scalar> fields.
// template<typename Dimension>
// void
// ConstantBoundary<Dimension>::
// enforceBoundary(Field<Dimension, std::vector<typename Dimension::Scalar> >& field) const {
//   resetValues(field, this->nodeIndices(), mVectorScalarValues);
// }

//------------------------------------------------------------------------------
// On problem startup we take our snapshot of the state.  We also then
// destroy the original internal nodes, as we will be replacing them
// with ghosts from now on.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConstantBoundary<Dimension>::initializeProblemStartup() {

  // Clear any existing data.
  mIntValues.clear();
  mScalarValues.clear();
  mVectorValues.clear();
  mTensorValues.clear();
  mSymTensorValues.clear();
  mThirdRankTensorValues.clear();
  mVectorScalarValues.clear();

  // Now take a snapshot of the Fields.
  const vector<int> nodeIDs = this->nodeIndices();
  storeFieldValues(*mNodeListPtr, nodeIDs, mIntValues);
  storeFieldValues(*mNodeListPtr, nodeIDs, mScalarValues);
  storeFieldValues(*mNodeListPtr, nodeIDs, mVectorValues);
  storeFieldValues(*mNodeListPtr, nodeIDs, mTensorValues);
  storeFieldValues(*mNodeListPtr, nodeIDs, mSymTensorValues);
  storeFieldValues(*mNodeListPtr, nodeIDs, mThirdRankTensorValues);
  storeFieldValues(*mNodeListPtr, nodeIDs, mVectorScalarValues);

  // Remove the origial internal nodes.
  mNodeListPtr->deleteNodes(nodeIDs);
  mNodeIDs.clear();

  // Turn the BC on.
  mActive = true;

  // Once we're done the boundary condition should be in a valid state.
  ENSURE(valid());
}

//------------------------------------------------------------------------------
// Test if the constant boundary is minimally valid.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
ConstantBoundary<Dimension>::valid() const {
  return true;
}

}
}

