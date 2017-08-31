//---------------------------------Spheral++----------------------------------//
// ConstantBoundary -- A boundary condition to enforce a constant 
// velocity on a given set of nodes.
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
#include "Boundary.hh"
#include "mapPositionThroughPlanes.hh"
#include "Field/FieldList.hh"
#include "Field/Field.hh"
#include "Field/FieldBase.hh"
#include "Utilities/planarReflectingOperator.hh"
#include "Utilities/DBC.hh"

#include "ConstantBoundary.hh"

#include "boost/lexical_cast.hpp"
#include "boost/algorithm/string.hpp"

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
using FileIOSpace::FileIO;

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
        CHECK2(*nodeItr >= 0 && *nodeItr < field.numElements(), *nodeItr << " " << field.numElements());
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
      // cerr << "Stored " << vals.size() << " values for " << key << endl;
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
    const std::vector<DataType>& newValues = itr->second;
    CHECK2(newValues.size() == nodeIDs.size(), key << " " << newValues.size() << " " << nodeIDs.size());
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
  mBoundaryCount(nodeList.numFields()),
  mNodeFlags("ConstantBoundaryNodeFlags" + boost::lexical_cast<std::string>(mBoundaryCount), nodeList, 0),
  mNumConstantNodes(nodeIDs.size()),
  mDenialPlane(denialPlane),
  mReflectOperator(planarReflectingOperator(denialPlane)),
  mActive(false),
  mScalarValues(),
  mVectorValues(),
  mTensorValues(),
  mSymTensorValues(),
  mThirdRankTensorValues(),
  mVectorScalarValues(),
  mVectorVectorValues(),
  mRestart(DataOutput::registerWithRestart(*this)) {

  // Store the ids of the nodes we're watching.
  for (vector<int>::const_iterator itr = nodeIDs.begin();
       itr < nodeIDs.end();
       ++itr) {
    REQUIRE(*itr >= 0.0 && *itr < nodeList.numInternalNodes());
    mNodeFlags[*itr] = 1;
  }

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

  if (mActive and &nodeList == mNodeListPtr) {
    typename Boundary<Dimension>::BoundaryNodes& boundaryNodes = this->accessBoundaryNodes(nodeList);
    vector<int>& cNodes = boundaryNodes.controlNodes;
    vector<int>& gNodes = boundaryNodes.ghostNodes;
    
    unsigned currentNumGhostNodes = nodeList.numGhostNodes();
    unsigned firstNewGhostNode = nodeList.numNodes();
    nodeList.numGhostNodes(currentNumGhostNodes + mNumConstantNodes);
    cNodes = vector<int>(mNumConstantNodes);
    gNodes = vector<int>(mNumConstantNodes);
    for (int i = 0; i != mNumConstantNodes; ++i) {
      const int j = firstNewGhostNode + i;
      mNodeFlags(j) = 1;
      cNodes[i] = j;
      gNodes[i] = j;
    }
    
    this->updateGhostNodes(nodeList);
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
  // const std::vector<int> ids = this->nodeIndices();
  // cerr << "Set ghosts for " << field.name() << " :";
  // for (std::vector<int>::const_iterator itr = ids.begin(); itr != ids.end(); ++itr) cerr << " " << field(*itr);
  // cerr << endl;
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

// Specialization for vector<Scalar> fields.
template<typename Dimension>
void
ConstantBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, std::vector<typename Dimension::Scalar>>& field) const {
  if (mActive) resetValues(field, this->nodeIndices(), mVectorScalarValues, false);
}

// Specialization for vector<Vector> fields.
template<typename Dimension>
void
ConstantBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, std::vector<typename Dimension::Vector>>& field) const {
  if (mActive) resetValues(field, this->nodeIndices(), mVectorVectorValues, false);
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

//------------------------------------------------------------------------------
// On problem startup we take our snapshot of the state.  We also then
// destroy the original internal nodes, as we will be replacing them
// with ghosts from now on.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConstantBoundary<Dimension>::initializeProblemStartup() {

  if (not mActive) {

    // Clear any existing data.
    mIntValues.clear();
    mScalarValues.clear();
    mVectorValues.clear();
    mTensorValues.clear();
    mSymTensorValues.clear();
    mThirdRankTensorValues.clear();
    mVectorScalarValues.clear();
    mVectorVectorValues.clear();

    // Now take a snapshot of the Fields.
    const vector<int> nodeIDs = this->nodeIndices();
    // cerr << "Node IDs: ";
    // std::copy(nodeIDs.begin(), nodeIDs.end(), std::ostream_iterator<int>(std::cerr, " "));
    // cerr << endl;
    storeFieldValues(*mNodeListPtr, nodeIDs, mIntValues);
    storeFieldValues(*mNodeListPtr, nodeIDs, mScalarValues);
    storeFieldValues(*mNodeListPtr, nodeIDs, mVectorValues);
    storeFieldValues(*mNodeListPtr, nodeIDs, mTensorValues);
    storeFieldValues(*mNodeListPtr, nodeIDs, mSymTensorValues);
    storeFieldValues(*mNodeListPtr, nodeIDs, mThirdRankTensorValues);
    storeFieldValues(*mNodeListPtr, nodeIDs, mVectorScalarValues);
    storeFieldValues(*mNodeListPtr, nodeIDs, mVectorVectorValues);

    // for (typename std::map<KeyType, std::vector<Scalar> >::const_iterator itr = mScalarValues.begin();
    //      itr != mScalarValues.end();
    //      ++itr) {
    //   cerr << " --> " << itr->first << " : ";
    //   std::copy(itr->second.begin(), itr->second.end(), std::ostream_iterator<double>(std::cerr, " "));
    //   cerr << endl;
    // }

    // Remove the origial internal nodes.
    mNodeListPtr->deleteNodes(nodeIDs);

    // Turn the BC on.
    mActive = true;
  }

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

//------------------------------------------------------------------------------
// Return a unique label for restart.
//------------------------------------------------------------------------------
template<typename Dimension>
std::string
ConstantBoundary<Dimension>::label() const {
  return "ConstantBoundary" + boost::lexical_cast<std::string>(mBoundaryCount);
}

//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConstantBoundary<Dimension>::
dumpState(FileIO& file, string pathName) const {
  file.write(mNumConstantNodes, pathName + "/numConstantNodes");
  file.write(mActive, pathName + "/active");
  file.write(mBoundaryCount, pathName + "/boundaryCount");

  vector<std::string> keys;
  for (const auto& p: mIntValues) {
    keys.push_back(p.first);
    file.write(p.second, pathName + "/IntValues/" + p.first);
  }
  file.write(keys, pathName + "/IntValues/keys");

  keys.clear();
  for (const auto& p: mScalarValues) {
    keys.push_back(p.first);
    file.write(p.second, pathName + "/ScalarValues/" + p.first);
  }
  file.write(keys, pathName + "/ScalarValues/keys");

  keys.clear();
  for (const auto& p: mVectorValues) {
    keys.push_back(p.first);
    file.write(p.second, pathName + "/VectorValues/" + p.first);
  }
  file.write(keys, pathName + "/VectorValues/keys");

  keys.clear();
  for (const auto& p: mTensorValues) {
    keys.push_back(p.first);
    file.write(p.second, pathName + "/TensorValues/" + p.first);
  }
  file.write(keys, pathName + "/TensorValues/keys");

  keys.clear();
  for (const auto& p: mSymTensorValues) {
    keys.push_back(p.first);
    file.write(p.second, pathName + "/SymTensorValues/" + p.first);
  }
  file.write(keys, pathName + "/SymTensorValues/keys");

  keys.clear();
  for (const auto& p: mThirdRankTensorValues) {
    keys.push_back(p.first);
    file.write(p.second, pathName + "/ThirdRankTensorValues/" + p.first);
  }
  file.write(keys, pathName + "/ThirdRankTensorValues/keys");

  keys.clear();
  for (const auto& p: mVectorScalarValues) {
    keys.push_back(p.first);
    file.write(p.second, pathName + "/VectorScalarValues/" + p.first);
  }
  file.write(keys, pathName + "/VectorScalarValues/keys");

  keys.clear();
  for (const auto& p: mVectorVectorValues) {
    keys.push_back(p.first);
    file.write(p.second, pathName + "/VectorVectorValues/" + p.first);
  }
  file.write(keys, pathName + "/VectorVectorValues/keys");
}

//------------------------------------------------------------------------------
// Read the state from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConstantBoundary<Dimension>::
restoreState(const FileIO& file, string pathName)  {
  file.read(mNumConstantNodes, pathName + "/numConstantNodes");
  file.read(mActive, pathName + "/active");
  file.read(mBoundaryCount, pathName + "/boundaryCount");

  vector<std::string> keys;
  file.read(keys, pathName + "/IntValues/keys");
  mIntValues.clear();
  for (const auto key: keys) {
    mIntValues[key] = std::vector<int>();
    file.read(mIntValues[key], pathName + "/IntValues/" + key);
  }

  keys.clear();
  file.read(keys, pathName + "/ScalarValues/keys");
  mScalarValues.clear();
  for (const auto key: keys) {
    mScalarValues[key] = std::vector<Scalar>();
    file.read(mScalarValues[key], pathName + "/ScalarValues/" + key);
  }

  keys.clear();
  file.read(keys, pathName + "/VectorValues/keys");
  mVectorValues.clear();
  for (const auto key: keys) {
    mVectorValues[key] = std::vector<Vector>();
    file.read(mVectorValues[key], pathName + "/VectorValues/" + key);
  }

  keys.clear();
  file.read(keys, pathName + "/TensorValues/keys");
  mTensorValues.clear();
  for (const auto key: keys) {
    mTensorValues[key] = std::vector<Tensor>();
    file.read(mTensorValues[key], pathName + "/TensorValues/" + key);
  }

  keys.clear();
  file.read(keys, pathName + "/SymTensorValues/keys");
  mSymTensorValues.clear();
  for (const auto key: keys) {
    mSymTensorValues[key] = std::vector<SymTensor>();
    file.read(mSymTensorValues[key], pathName + "/SymTensorValues/" + key);
  }

  keys.clear();
  file.read(keys, pathName + "/ThirdRankTensorValues/keys");
  mThirdRankTensorValues.clear();
  for (const auto key: keys) {
    mThirdRankTensorValues[key] = std::vector<ThirdRankTensor>();
    file.read(mThirdRankTensorValues[key], pathName + "/ThirdRankTensorValues/" + key);
  }

  keys.clear();
  file.read(keys, pathName + "/VectorScalarValues/keys");
  mVectorScalarValues.clear();
  for (const auto key: keys) {
    mVectorScalarValues[key] = std::vector<std::vector<Scalar> >();
    file.read(mVectorScalarValues[key], pathName + "/VectorScalarValues/" + key);
  }

  keys.clear();
  file.read(keys, pathName + "/VectorVectorValues/keys");
  mVectorVectorValues.clear();
  for (const auto key: keys) {
    mVectorVectorValues[key] = std::vector<std::vector<Vector> >();
    file.read(mVectorVectorValues[key], pathName + "/VectorSVectorValues/" + key);
  }
}

}
}

