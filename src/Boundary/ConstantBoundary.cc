//---------------------------------Spheral++----------------------------------//
// ConstantBoundary -- A boundary condition to enforce a constant 
// velocity on a given set of nodes.
//----------------------------------------------------------------------------//
#include "ConstantBoundary.hh"
#include "Boundary.hh"
#include "Field/FieldList.hh"
#include "Field/Field.hh"
#include "Field/FieldBase.hh"

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
      const std::string key = StateBase<Dimension>::key(**fieldItr);
      CHECK(values.find(key) == values.end());
      values[key] = vals;
    }
  }
}
  
//------------------------------------------------------------------------------
// Set the ghost values in the given field using the given map.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
void
setGhostValues(FieldSpace::Field<Dimension, DataType>& field,
               const std::vector<int>& ghostNodes,
               const std::map<std::string, std::vector<DataType> >& values) {

  const NodeSpace::NodeList<Dimension>& nodeList = field.nodeList();

  // Find this Field in the set of stored values.
  const std::string key = StateBase<Dimension>::key(field);
  typename std::map<std::string, std::vector<DataType> >::const_iterator itr = values.find(key);

  // Now set the ghost values.
  if (itr != values.end()) {
    const std::vector<DataType>& ghostValues = itr->second;;
    CHECK(ghostValues.size() == ghostNodes.size());
    for (int i = 0; i < ghostNodes.size(); ++i) {
      CHECK(ghostNodes[i] >= nodeList.firstGhostNode() &&
            ghostNodes[i] < nodeList.numNodes());
      field(ghostNodes[i]) = ghostValues[i];
    }
  }

}

}

//------------------------------------------------------------------------------
// Construct with the given set of nodes.
//------------------------------------------------------------------------------
template<typename Dimension>
ConstantBoundary<Dimension>::
ConstantBoundary(const NodeList<Dimension>& nodeList,
                 const vector<int>& nodeIDs):
  Boundary<Dimension>(),
  mNodeListPtr(&nodeList),
  mNodeIDs(nodeIDs),
  mNumConstantNodes(nodeIDs.size()),
  mScalarValues(),
  mVectorValues(),
  mTensorValues(),
  mSymTensorValues(),
  mThirdRankTensorValues(),
  mVectorScalarValues() {

  // Make an initial copy of the ghost state.
  this->initializeProblemStartup();

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
    setGhostValues(positions, this->ghostNodes(nodeList), mVectorValues);

    // Set the ghost H fields.
    Field<Dimension, SymTensor>& Hfield = nodeList.Hfield();
    setGhostValues(Hfield, this->ghostNodes(nodeList), mSymTensorValues);

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
  setGhostValues(field, this->ghostNodes(field.nodeList()), mIntValues);
}

// Specialization for scalar fields.
template<typename Dimension>
void
ConstantBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::Scalar>& field) const {
  cerr << "--> " << field.name() << " ";
  const std::string key = StateBase<Dimension>::key(field);
  typename std::map<std::string, std::vector<Scalar> >::const_iterator itr = mScalarValues.find(key);
  cerr << (itr != mScalarValues.end()) << " : ";
  // if (itr != mScalarValues.end()) {
  //   std::copy(itr->second.begin(), itr->second.end(), std::ostream_iterator<Scalar>(std::cerr, " "));
  // }
  std::cerr << endl;
  setGhostValues(field, this->ghostNodes(field.nodeList()), mScalarValues);
}

// Specialization for Vector fields.
template<typename Dimension>
void
ConstantBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::Vector>& field) const {
  setGhostValues(field, this->ghostNodes(field.nodeList()), mVectorValues);
}

// Specialization for Tensor fields.
template<typename Dimension>
void
ConstantBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::Tensor>& field) const {
  setGhostValues(field, this->ghostNodes(field.nodeList()), mTensorValues);
}

// Specialization for symmetric tensors.
template<typename Dimension>
void
ConstantBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::SymTensor>& field) const {
  setGhostValues(field, this->ghostNodes(field.nodeList()), mSymTensorValues);
}

// Specialization for third rank tensors.
template<typename Dimension>
void
ConstantBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::ThirdRankTensor>& field) const {
  setGhostValues(field, this->ghostNodes(field.nodeList()), mThirdRankTensorValues);
}

// Specialization for vector<scalar> fields.
template<typename Dimension>
void
ConstantBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, std::vector<typename Dimension::Scalar> >& field) const {
  setGhostValues(field, this->ghostNodes(field.nodeList()), mVectorScalarValues);
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
// On problem startup we take our snapshot of the state.
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
  storeFieldValues(*mNodeListPtr, mNodeIDs, mIntValues);
  storeFieldValues(*mNodeListPtr, mNodeIDs, mScalarValues);
  storeFieldValues(*mNodeListPtr, mNodeIDs, mVectorValues);
  storeFieldValues(*mNodeListPtr, mNodeIDs, mTensorValues);
  storeFieldValues(*mNodeListPtr, mNodeIDs, mSymTensorValues);
  storeFieldValues(*mNodeListPtr, mNodeIDs, mThirdRankTensorValues);
  storeFieldValues(*mNodeListPtr, mNodeIDs, mVectorScalarValues);

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

