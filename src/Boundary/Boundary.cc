//---------------------------------Spheral++----------------------------------//
// Boundary -- Abstract base class for the boundary condition classes.
//
// Created by JMO, Wed Feb 16 21:01:02 PST 2000
//----------------------------------------------------------------------------//

#include <algorithm>

#include "Boundary.hh"
#include "DataBase/DataBase.hh"
#include "NodeList/NodeList.hh"
#include "Field/FieldList.hh"
#include "Utilities/removeElements.hh"

#include "Utilities/DBC.hh"

using std::vector;
using std::map;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {

//------------------------------------------------------------------------------
// Empty constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
Boundary<Dimension>::Boundary():
  mBoundaryNodes() {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
Boundary<Dimension>::~Boundary() {
}

//------------------------------------------------------------------------------
// Allow read access to the map of NodeList->BoundaryNodes.
//------------------------------------------------------------------------------
template<typename Dimension>
const map<NodeList<Dimension>*, typename Boundary<Dimension>::BoundaryNodes>&
Boundary<Dimension>::boundaryNodeMap() const {
  return mBoundaryNodes;
}

//------------------------------------------------------------------------------
// Return the list of control nodes for the given NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
const vector<int>&
Boundary<Dimension>::controlNodes(const NodeList<Dimension>& nodeList) const {
  typename map<NodeList<Dimension>*, BoundaryNodes>::const_iterator
    itr = mBoundaryNodes.find(const_cast<NodeList<Dimension>*>(&nodeList));
  if (itr != mBoundaryNodes.end()) {
    return itr->second.controlNodes;
  } else {
    VERIFY2(false, "Boundary::controlNodes: no entry for NodeList: " << nodeList.name());
  }
}

//------------------------------------------------------------------------------
// Return the list of ghost nodes for the given NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
const vector<int>&
Boundary<Dimension>::ghostNodes(const NodeList<Dimension>& nodeList) const {
  typename map<NodeList<Dimension>*, BoundaryNodes>::const_iterator 
    itr = mBoundaryNodes.find(const_cast<NodeList<Dimension>*>(&nodeList));
  if (itr != mBoundaryNodes.end()) {
    return itr->second.ghostNodes;
  } else {
    std::cerr << "Number of known NodeLists: " << mBoundaryNodes.size() << std::endl;
    VERIFY2(false, "Boundary::ghostNodes: no entry for NodeList: " << nodeList.name());
  }
}

//------------------------------------------------------------------------------
// Return the list of violation nodes for the given NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
const vector<int>&
Boundary<Dimension>::violationNodes(const NodeList<Dimension>& nodeList) const {
  typename map<NodeList<Dimension>*, BoundaryNodes>::const_iterator 
    itr = mBoundaryNodes.find(const_cast<NodeList<Dimension>*>(&nodeList));
  if (itr != mBoundaryNodes.end()) {
    return itr->second.violationNodes;
  } else {
    VERIFY2(false, "Boundary::violationNodes: no entry for NodeList: " << nodeList.name());
  }
}

//------------------------------------------------------------------------------
// Return begin/end iterators for the control nodes of the given NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
typename vector<int>::const_iterator
Boundary<Dimension>::controlBegin(const NodeList<Dimension>& nodeList) const {
  return controlNodes(nodeList).begin();
}

template<typename Dimension>
typename vector<int>::const_iterator
Boundary<Dimension>::controlEnd(const NodeList<Dimension>& nodeList) const {
  return controlNodes(nodeList).end();
}

//------------------------------------------------------------------------------
// Return begin/end iterators for the ghost nodes of the given NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
typename vector<int>::const_iterator
Boundary<Dimension>::ghostBegin(const NodeList<Dimension>& nodeList) const {
  return ghostNodes(nodeList).begin();
}

template<typename Dimension>
typename vector<int>::const_iterator
Boundary<Dimension>::ghostEnd(const NodeList<Dimension>& nodeList) const {
  return ghostNodes(nodeList).end();
}

//------------------------------------------------------------------------------
// Return begin/end iterators for the violation nodes of the given NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
typename vector<int>::const_iterator
Boundary<Dimension>::violationBegin(const NodeList<Dimension>& nodeList) const {
  return violationNodes(nodeList).begin();
}

template<typename Dimension>
typename vector<int>::const_iterator
Boundary<Dimension>::violationEnd(const NodeList<Dimension>& nodeList) const {
  return violationNodes(nodeList).end();
}

//------------------------------------------------------------------------------
// Default method to set ghost nodes for all NodeLists in DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Boundary<Dimension>::setAllGhostNodes(DataBase<Dimension>& dataBase) {
  for (typename DataBase<Dimension>::NodeListIterator 
	 nodeListItr = dataBase.nodeListBegin();
       nodeListItr < dataBase.nodeListEnd();
       ++nodeListItr) {
    setGhostNodes(**nodeListItr);
  }
}

//------------------------------------------------------------------------------
// Default method to set violation nodes for all NodeLists in DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Boundary<Dimension>::setAllViolationNodes(DataBase<Dimension>& dataBase) {
  for (typename DataBase<Dimension>::NodeListIterator 
	 nodeListItr = dataBase.nodeListBegin();
       nodeListItr < dataBase.nodeListEnd();
       ++nodeListItr) {
    setViolationNodes(**nodeListItr);
  }

  BEGIN_CONTRACT_SCOPE
  // Make sure the boundary knows about each nodelist in the database.
  for (typename DataBase<Dimension>::ConstNodeListIterator 
       i = dataBase.nodeListBegin(); i != dataBase.nodeListEnd(); ++i)
  {
    ENSURE(mBoundaryNodes.find(*i) != mBoundaryNodes.end());
  } // end for
  END_CONTRACT_SCOPE
}

//------------------------------------------------------------------------------
// Cull out inactive ghost nodes based on a FieldList of flags.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Boundary<Dimension>::cullGhostNodes(const FieldList<Dimension, int>& flagSet,
                                    FieldList<Dimension, int>& old2newIndexMap,
                                    vector<int>& numNodesRemoved) {

  typedef NodeListRegistrar<Dimension> Registrar;
  Registrar& registrar = Registrar::instance();
  REQUIRE(numNodesRemoved.size() == registrar.numNodeLists());

  // Walk the NodeLists.
  size_t nodeListi = 0;
  for (typename Registrar::const_iterator itr = registrar.begin();
       itr != registrar.end();
       ++itr, ++nodeListi) {
    const NodeList<Dimension>* nodeListPtr = *itr;

    // Does the Boundary have entries for this NodeList?
    if (haveNodeList(*nodeListPtr)) {
      BoundaryNodes& boundaryNodes = accessBoundaryNodes(const_cast<NodeList<Dimension>&>(*nodeListPtr));
      CHECK(boundaryNodes.ghostNodes.size() == boundaryNodes.controlNodes.size());
      if (boundaryNodes.ghostNodes.size() > 0) {
        const size_t myOldFirstGhostNode = boundaryNodes.ghostNodes[0];
        const size_t myNewFirstGhostNode = myOldFirstGhostNode - numNodesRemoved[nodeListi];

        // Grab the flags for this NodeList.
        CHECK(flagSet.haveNodeList(*nodeListPtr));
        const Field<Dimension, int>& flags = *(flagSet[nodeListi]);

        // Patch up the ghost and control node indices.
        vector<int> newGhostNodes, newControlNodes;
        size_t newGhostIndex = myNewFirstGhostNode;
        for (size_t k = 0; k != boundaryNodes.ghostNodes.size(); ++k) {
          if (flags(boundaryNodes.ghostNodes[k]) == 1) {
            newGhostNodes.push_back(newGhostIndex);
            old2newIndexMap(nodeListi, boundaryNodes.ghostNodes[k]) = newGhostIndex;
            ++newGhostIndex;
            // CHECK(boundaryNodes.controlNodes[k] < myOldFirstGhostNode);  // <-- Not true for the ConstantBoundary
            newControlNodes.push_back(old2newIndexMap(nodeListi, boundaryNodes.controlNodes[k]));
            // CHECK(newControlNodes.back() < myNewFirstGhostNode);  // <-- Not true for the ConstantBoundary
          } else {
            ++numNodesRemoved[nodeListi];
          }
        }

        // Update the ghost & control nodes, and the result indicating how many of our ghost nodes 
        // were removed.
        boundaryNodes.ghostNodes = newGhostNodes;
        boundaryNodes.controlNodes = newControlNodes;
        CHECK(boundaryNodes.ghostNodes.size() == boundaryNodes.controlNodes.size());
      }
    }
  }
}

//------------------------------------------------------------------------------
// Return a non-const reference to the full map of BoundaryNodes.
//------------------------------------------------------------------------------
template<typename Dimension>
map<NodeList<Dimension>*, typename Boundary<Dimension>::BoundaryNodes>&
Boundary<Dimension>::accessBoundaryNodes() {
  return mBoundaryNodes;
}

//------------------------------------------------------------------------------
// Return a non-const reference to the current list of control nodes.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Boundary<Dimension>::BoundaryNodes&
Boundary<Dimension>::accessBoundaryNodes(NodeList<Dimension>& nodeList) {
  typename map<NodeList<Dimension>*, BoundaryNodes>::const_iterator 
    itr = mBoundaryNodes.find(&nodeList);
  if (itr != mBoundaryNodes.end()) {
    return mBoundaryNodes[&nodeList];
  } else {
    VERIFY2(false, "Boundary::accessBoundaryNodes: no entry for NodeList: " << nodeList.name());
  }
}

//------------------------------------------------------------------------------
// Add the given NodeList to the current set of BoundaryNodes.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Boundary<Dimension>::addNodeList(NodeList<Dimension>& nodeList) {
  if (mBoundaryNodes.find(&nodeList) == mBoundaryNodes.end()) {
    mBoundaryNodes[&nodeList] = BoundaryNodes();
  }
}

//------------------------------------------------------------------------------
// Default for vector<Scalar> fields, just perform a copy.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Boundary<Dimension>::
applyGhostBoundary(Field<Dimension, std::vector<typename Dimension::Scalar>>& field) const {
  const auto& nodeList = field.nodeList();
  CHECK(this->controlNodes(nodeList).size() == this->ghostNodes(nodeList).size());
  auto controlItr = this->controlBegin(nodeList);
  auto ghostItr = this->ghostBegin(nodeList);
  for (; controlItr < this->controlEnd(nodeList); ++controlItr, ++ghostItr) {
    CHECK(ghostItr < this->ghostEnd(nodeList));
    CHECK(*controlItr >= 0 && *controlItr < nodeList.numNodes());
    CHECK(*ghostItr >= nodeList.firstGhostNode() && *ghostItr < nodeList.numNodes());
    field(*ghostItr) = field(*controlItr);
  }
}

//------------------------------------------------------------------------------
// Default for vector<Vector> fields, just perform a copy
//------------------------------------------------------------------------------
template<typename Dimension>
void
Boundary<Dimension>::
applyGhostBoundary(Field<Dimension, std::vector<typename Dimension::Vector> >& field) const {
  const auto& nodeList = field.nodeList();
  CHECK(this->controlNodes(nodeList).size() == this->ghostNodes(nodeList).size());
  auto controlItr = this->controlBegin(nodeList);
  auto ghostItr = this->ghostBegin(nodeList);
  for (; controlItr < this->controlEnd(nodeList); ++controlItr, ++ghostItr) {
    CHECK(ghostItr < this->ghostEnd(nodeList));
    CHECK(*controlItr >= 0 && *controlItr < nodeList.numNodes());
    CHECK(*ghostItr >= nodeList.firstGhostNode() && *ghostItr < nodeList.numNodes());
    field(*ghostItr) = field(*controlItr);
  }
}

//------------------------------------------------------------------------------
// Clear out any NodeList information that is currently present.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Boundary<Dimension>::reset(const DataBase<Dimension>& dataBase) {
  // Clear our own internal data.
  mBoundaryNodes = std::map<NodeList<Dimension>*, BoundaryNodes>();
}

//------------------------------------------------------------------------------
// Report the number of ghost nodes.
//------------------------------------------------------------------------------
template<typename Dimension>
int
Boundary<Dimension>::numGhostNodes() const {
  int result = 0;
  for (typename map<NodeList<Dimension>*, BoundaryNodes>::const_iterator itr = mBoundaryNodes.begin();
       itr != mBoundaryNodes.end();
       ++itr) {
    result += itr->second.ghostNodes.size();
  }
  return result;
}

//------------------------------------------------------------------------------
// Default clipping opertion to a no-op.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Boundary<Dimension>::
clip(typename Dimension::Vector& xmin, typename Dimension::Vector& xmax) const {
}

}
