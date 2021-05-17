//---------------------------------Spheral++----------------------------------//
// HostCodeBoundary -- A boundary condition for user-specified ghost nodes
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
#include "Boundary.hh"
#include "Field/FieldList.hh"
#include "Field/Field.hh"
#include "Field/FieldBase.hh"
#include "Utilities/DBC.hh"

#include "Boundary/HostCodeBoundary.hh"
#include "Boundary/ConstantBoundaryUtilities.hh"

#include "boost/algorithm/string.hpp"

using std::vector;
using std::map;
using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {

//------------------------------------------------------------------------------
// Construct with the given set of nodes.
//------------------------------------------------------------------------------
template<typename Dimension>
HostCodeBoundary<Dimension>::
HostCodeBoundary(DataBase<Dimension>& dataBase,
                 NodeList<Dimension>& nodeList,
                 const vector<int>& nodeIDs):
  Boundary<Dimension>(),
  mDataBase(dataBase),
  mNodeListPtr(&nodeList),
  mBoundaryCount(nodeList.numFields()),
  mNodeFlags("HostCodeBoundaryNodeFlags" + std::to_string(mBoundaryCount), nodeList, 0),
  mNumConstantNodes(nodeIDs.size()),
  mActive(false),
  mBufferedValues(),
  mRestart(registerWithRestart(*this)),
  mRedistribution(registerWithRedistribution(*this,
                                             &HostCodeBoundary<Dimension>::notifyBeforeRedistribution,
                                             &HostCodeBoundary<Dimension>::notifyAfterRedistribution)) {

  // Store the ids of the nodes we're watching.
  for (auto itr = nodeIDs.begin(); itr < nodeIDs.end(); ++itr) {
    REQUIRE(*itr >= 0.0 && *itr < (int)nodeList.numInternalNodes());
    mNodeFlags[*itr] = 1;
  }
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
HostCodeBoundary<Dimension>::~HostCodeBoundary() {
}

//------------------------------------------------------------------------------
// setGhostNodes
//------------------------------------------------------------------------------
template<typename Dimension>
void
HostCodeBoundary<Dimension>::
setGhostNodes(NodeList<Dimension>& nodeList) {
  this->addNodeList(nodeList);

  if (mActive and &nodeList == mNodeListPtr) {
    auto& boundaryNodes = this->accessBoundaryNodes(nodeList);
    auto& cNodes = boundaryNodes.controlNodes;
    auto& gNodes = boundaryNodes.ghostNodes;
    
    auto currentNumGhostNodes = nodeList.numGhostNodes();
    auto firstNewGhostNode = nodeList.numNodes();
    nodeList.numGhostNodes(currentNumGhostNodes + mNumConstantNodes);
    cNodes.resize(mNumConstantNodes);
    gNodes.resize(mNumConstantNodes);
    for (auto i = 0; i < mNumConstantNodes; ++i) {
      const auto j = firstNewGhostNode + i;
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
HostCodeBoundary<Dimension>::
updateGhostNodes(NodeList<Dimension>& nodeList) {
  if (mActive and &nodeList == mNodeListPtr) {
    this->applyGhostBoundary(nodeList.positions());
    this->applyGhostBoundary(nodeList.Hfield());
  }
}

//------------------------------------------------------------------------------
// Apply the boundary condition to the ghost nodes
//------------------------------------------------------------------------------
template<typename Dimension>
void
HostCodeBoundary<Dimension>::
applyGhostBoundary(FieldBase<Dimension>& field) const {
  if (mActive) resetValues(field, this->nodeIndices(), mBufferedValues, false);
}

//------------------------------------------------------------------------------
// setViolationNodes
//------------------------------------------------------------------------------
template<typename Dimension>
void
HostCodeBoundary<Dimension>::
setViolationNodes(NodeList<Dimension>& nodeList) {
}

//------------------------------------------------------------------------------
// updateViolationNodes, correcting positions and H's.
//------------------------------------------------------------------------------
template<typename Dimension>
void
HostCodeBoundary<Dimension>::
updateViolationNodes(NodeList<Dimension>& nodeList) {
}

//------------------------------------------------------------------------------
// enforceBoundary
//------------------------------------------------------------------------------
template<typename Dimension>
void
HostCodeBoundary<Dimension>::
enforceBoundary(FieldBase<Dimension>&) const {
  // resetValues(field, this->nodeIndices(), mIntValues);
}

//------------------------------------------------------------------------------
// On problem startup we take our snapshot of the state.  We also then
// destroy the original internal nodes, as we will be replacing them
// with ghosts from now on.
//------------------------------------------------------------------------------
template<typename Dimension>
void
HostCodeBoundary<Dimension>::initializeProblemStartup(const bool final) {

  // Clear any existing data.
  mBufferedValues.clear();

  // Now take a snapshot of the Fields.
  const auto nodeIDs = this->nodeIndices();
  // cerr << "Node IDs: ";
  // std::copy(nodeIDs.begin(), nodeIDs.end(), std::ostream_iterator<int>(std::cerr, " "));
  // cerr << endl;
  storeFieldValues(*mNodeListPtr, nodeIDs, mBufferedValues);

  // If we're in cylindrical symmetry (RZ) we need to convert the mass to mass/(2*pi*r).
  if (mDataBase.isRZ) {
    const auto& pos = mNodeListPtr->positions();
    const auto& mass = mNodeListPtr->mass();
    std::vector<char> buffer;
    for (auto i: nodeIDs) {
      const auto circi = 2.0*M_PI*abs(pos(i).y());
      CHECK(circi > 0.0);
      const auto mi = mass(i)/circi;
      packElement(mi, buffer);
    }
    const auto mkey = StateBase<Dimension>::key(mass);
    CHECK(mBufferedValues.find(mkey) != mBufferedValues.end());
    mBufferedValues[mkey] = buffer;
  }

  // Remove the origial internal nodes.
  //if (final) {
  //  mNodeListPtr->deleteNodes(nodeIDs);
  //  mActive = true;
  //}
}

//------------------------------------------------------------------------------
// updateNodeIDs
//------------------------------------------------------------------------------
template<typename Dimension>
void
HostCodeBoundary<Dimension>::
updateNodeIDs(const vector<int>& nodeIDs) {
  mNumConstantNodes = nodeIDs.size();
  // Zero out the previous flags
  for (int i = 0; i < mNodeListPtr->numNodes(); ++i) {
    mNodeFlags[i] = 0;
  }
  // Store the ids of the nodes we're watching.
  for (auto itr = nodeIDs.begin(); itr < nodeIDs.end(); ++itr) {
    mNodeFlags[*itr] = 1;
  }
  // call initializeProblemStartup
  this->initializeProblemStartup(true);
}

//------------------------------------------------------------------------------
// Return the set of node IDs we're controlling.
//------------------------------------------------------------------------------
template<typename Dimension>
std::vector<int>
HostCodeBoundary<Dimension>::
nodeIndices() const {
  std::vector<int> result;
  if (mActive) {
    for (const auto& bnitr: this->accessBoundaryNodes()) {
      const auto& ghosts = bnitr.second.ghostNodes;
      std::copy(ghosts.begin(), ghosts.end(), std::back_inserter(result));
    }
  } else {
    for (auto i = 0u; i != mNodeListPtr->numNodes(); ++i) {
      if (mNodeFlags(i) == 1) result.push_back(i);
    }
  }
  return result;
}

//------------------------------------------------------------------------------
// Return a unique label for restart.
//------------------------------------------------------------------------------
template<typename Dimension>
std::string
HostCodeBoundary<Dimension>::label() const {
  return "HostCodeBoundary" + std::to_string(mBoundaryCount);
}

//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
HostCodeBoundary<Dimension>::
dumpState(FileIO& file, const string& pathName) const {
  file.write(mNumConstantNodes, pathName + "/numConstantNodes");
  file.write(mActive, pathName + "/active");
  file.write(mBoundaryCount, pathName + "/boundaryCount");

  vector<std::string> keys;
  for (const auto& p: mBufferedValues) {
    keys.push_back(p.first);
    std::string val(p.second.begin(), p.second.end());
    file.write(val, pathName + "/BufferedValues/" + p.first);
  }
  file.write(keys, pathName + "/keys");
}

//------------------------------------------------------------------------------
// Read the state from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
HostCodeBoundary<Dimension>::
restoreState(const FileIO& file, const string& pathName)  {
  file.read(mNumConstantNodes, pathName + "/numConstantNodes");
  file.read(mActive, pathName + "/active");
  file.read(mBoundaryCount, pathName + "/boundaryCount");

  vector<std::string> keys;
  file.read(keys, pathName + "/keys");
  mBufferedValues.clear();
  for (const auto key: keys) {
    std::string val;
    file.read(val, pathName + "/BufferedValues/" + key);
    mBufferedValues[key] = vector<char>(val.begin(), val.end());
  }
}

//------------------------------------------------------------------------------
// Redistribution methods
//------------------------------------------------------------------------------
template<typename Dimension>
void
HostCodeBoundary<Dimension>::
notifyBeforeRedistribution() {
  VERIFY2(false, "HostCodeBoundary ERROR: node redistribution not allowed with constant boundaries.");
}

template<typename Dimension>
void
HostCodeBoundary<Dimension>::
notifyAfterRedistribution() {
  VERIFY2(false, "HostCodeBoundary ERROR: node redistribution not allowed with constant boundaries.");
}

}
