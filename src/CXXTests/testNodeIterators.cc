//------------------------------------------------------------------------------
// testNodeIterators
// A collection of test functions for the node iterators.
//
// Created by JMO, Thu Mar 25 14:21:41 2004
//------------------------------------------------------------------------------
#include "testNodeIterators.hh"

#include "Field/NodeIterators.hh"
#include "DataBase/DataBase.hh"
#include "NodeList/NodeList.hh"
#include "Utilities/DBC.hh"

#include "Geometry/Dimension.hh"

#include <algorithm>
#include <vector>
#include <map>
#include <string>

namespace Spheral {

using std::vector;
using std::map;
using std::string;
  
//------------------------------------------------------------------------------
// Workhorse method to test that a given pair of NodeIterators walks the given 
// set of node IDs.
//------------------------------------------------------------------------------
template<typename Dimension, typename IteratorType>
inline
string
testNodeIteratorPair(IteratorType beginIterator,
                     IteratorType endIterator,
                     const map<const NodeList<Dimension>*, vector<int> >& controlIDs) {

  // Set up a list of flags for each control ID, to make sure we hit each one.
  map<const NodeList<Dimension>*, vector<bool> > flags;
  int numCache = 0;
  for (typename map<const NodeList<Dimension>*, vector<int> >::const_iterator itr = controlIDs.begin();
       itr != controlIDs.end();
       ++itr) {
    flags[itr->first] = vector<bool>(itr->second.size(), false);
    numCache += itr->second.size();
  }
  CHECK(flags.size() == controlIDs.size());

  // Now go over the iterators and flag each node we hit.
  int lastCacheID = -1;
  for (IteratorType itr = beginIterator; itr != endIterator; ++itr) {
//     cerr << itr.nodeID() << " "
//          << itr.fieldID() << " "
//          << itr.cacheID() << " : "
//          << endIterator.nodeID() << " "
//          << endIterator.fieldID() << " "
//          << endIterator.cacheID() << " : "
//          << (itr == endIterator) << endl;

//     // Check the cache information.
//     if (itr.numCache() != numCache) {
//       return "Incorrect number of cache IDs in iterator.";
//     }
//     if (itr.cacheID() != lastCacheID + 1) {
//       return "Incorrect cache ID.";
//     }
//     ++lastCacheID;

    // Check that this NodeList is in the expected set.
    const NodeList<Dimension>* nodeListPtr = itr.nodeListPtr();
    const int nodeID = itr.nodeID();
    if (controlIDs.find(nodeListPtr) == controlIDs.end()) {
      return "Iterator contains invalid NodeList";
    }

    // Check that the node ID is in the expected set for this NodeList.
    const vector<int>& cIDs = controlIDs.find(nodeListPtr)->second;
    vector<int>::const_iterator IDItr = find(cIDs.begin(), cIDs.end(), nodeID);
    if (IDItr == cIDs.end()) {
      return "Iterator contains invalid node ID for valid NodeList.";
    }

    // Check whether this node ID has already been hit.
    const int i = distance(cIDs.begin(), IDItr);
    CHECK(flags.find(nodeListPtr) != flags.end());
    CHECK(i < flags[nodeListPtr].size());
    if (flags[nodeListPtr][i] == true) {
      return "Iterator contains redundant node ID for valid Node and NodeList";
    }

    // Flag the node ID as hit.
    flags[nodeListPtr][i] = true;

  }

  // Check that we hit all node IDs in the control set.
  for (typename map<const NodeList<Dimension>*, vector<bool> >::const_iterator itr = flags.begin();
       itr != flags.end();
       ++itr) {
    if (find(itr->second.begin(), itr->second.end(), false) != itr->second.end()) {
      return "At least one of the expected control nodes was not walked by the iterator.";
    }
  }

  // Check that the end iterator has the expected state.
  if (endIterator.nodeID() != 0) {
    return "End iterator does not have expected nodeID == 0";
  }
//   if (endIterator.fieldID() != std::distance(beginIterator.nodeListIterator(), 
//                                              endIterator.nodeListIterator())) {
//     return "End iterator does not have expected fieldID";
//   }
  if (endIterator.nodeListPtr() != 0) {
    return "End iterator does not have expected nodeListPtr == 0";
  }

  // If we made it to here, everything checks out.
  return "OK";
}

//------------------------------------------------------------------------------
// Test global AllNodeIterators.
//------------------------------------------------------------------------------
template<typename Dimension>
std::string
testGlobalAllNodeIterators(const DataBase<Dimension>& dataBase) {

  // Build the set of node IDs we expect the iterators to span.
  std::map<const NodeList<Dimension>*, std::vector<int> > controlIDs;
  for (typename DataBase<Dimension>::ConstNodeListIterator nodeListItr = dataBase.nodeListBegin();
       nodeListItr != dataBase.nodeListEnd();
       ++nodeListItr) {
    controlIDs[*nodeListItr] = std::vector<int>();
    controlIDs[*nodeListItr].reserve((*nodeListItr)->numNodes());
    for (int i = 0; i != (*nodeListItr)->numNodes(); ++i) controlIDs[*nodeListItr].push_back(i);
    CHECK(controlIDs[*nodeListItr].size() == (*nodeListItr)->numNodes());
  }
  CHECK(controlIDs.size() == dataBase.numNodeLists());

  // Perform the check.
  return testNodeIteratorPair(dataBase.nodeBegin(),
                              dataBase.nodeEnd(),
                              controlIDs);
}

//------------------------------------------------------------------------------
// Test global InternalNodeIterators.
//------------------------------------------------------------------------------
template<typename Dimension>
std::string
testGlobalInternalNodeIterators(const DataBase<Dimension>& dataBase) {

  // Build the set of node IDs we expect the iterators to span.
  std::map<const NodeList<Dimension>*, std::vector<int> > controlIDs;
  for (typename DataBase<Dimension>::ConstNodeListIterator nodeListItr = dataBase.nodeListBegin();
       nodeListItr != dataBase.nodeListEnd();
       ++nodeListItr) {
    controlIDs[*nodeListItr] = std::vector<int>();
    controlIDs[*nodeListItr].reserve((*nodeListItr)->numInternalNodes());
    for (int i = 0; i != (*nodeListItr)->numInternalNodes(); ++i) controlIDs[*nodeListItr].push_back(i);
    CHECK(controlIDs[*nodeListItr].size() == (*nodeListItr)->numInternalNodes());
  }
  CHECK(controlIDs.size() == dataBase.numNodeLists());

  // Perform the check.
  return testNodeIteratorPair(dataBase.internalNodeBegin(),
                              dataBase.internalNodeEnd(),
                              controlIDs);
}

//------------------------------------------------------------------------------
// Test global GhostNodeIterators.
//------------------------------------------------------------------------------
template<typename Dimension>
std::string
testGlobalGhostNodeIterators(const DataBase<Dimension>& dataBase) {

  // Build the set of node IDs we expect the iterators to span.
  std::map<const NodeList<Dimension>*, std::vector<int> > controlIDs;
  for (typename DataBase<Dimension>::ConstNodeListIterator nodeListItr = dataBase.nodeListBegin();
       nodeListItr != dataBase.nodeListEnd();
       ++nodeListItr) {
    controlIDs[*nodeListItr] = std::vector<int>();
    controlIDs[*nodeListItr].reserve((*nodeListItr)->numGhostNodes());
    for (int i = (*nodeListItr)->firstGhostNode(); i != (*nodeListItr)->numNodes(); ++i) controlIDs[*nodeListItr].push_back(i);
    CHECK(controlIDs[*nodeListItr].size() == (*nodeListItr)->numGhostNodes());
  }
  CHECK(controlIDs.size() == dataBase.numNodeLists());

  // Perform the check.
  return testNodeIteratorPair(dataBase.ghostNodeBegin(),
                              dataBase.ghostNodeEnd(),
                              controlIDs);
}

//------------------------------------------------------------------------------
// Test global MasterNodeIterators.
//------------------------------------------------------------------------------
template<typename Dimension>
std::string
testGlobalMasterNodeIterators(const DataBase<Dimension>& dataBase) {

  // Build the set of node IDs we expect the iterators to span.
  std::map<const NodeList<Dimension>*, std::vector<int> > controlIDs;
  for (typename DataBase<Dimension>::ConstNodeListIterator nodeListItr = dataBase.nodeListBegin();
       nodeListItr != dataBase.nodeListEnd();
       ++nodeListItr) {
    controlIDs[*nodeListItr] = (*nodeListItr)->neighbor().masterList();
    CHECK(controlIDs[*nodeListItr].size() == (*nodeListItr)->neighbor().numMaster());
  }
  CHECK(controlIDs.size() == dataBase.numNodeLists());

  // Perform the check.
  return testNodeIteratorPair(dataBase.masterNodeBegin(),
                              dataBase.masterNodeEnd(),
                              controlIDs);
}

//------------------------------------------------------------------------------
// Test global CoarseNodeIterators.
//------------------------------------------------------------------------------
template<typename Dimension>
std::string
testGlobalCoarseNodeIterators(const DataBase<Dimension>& dataBase) {

  // Build the set of node IDs we expect the iterators to span.
  std::map<const NodeList<Dimension>*, std::vector<int> > controlIDs;
  for (typename DataBase<Dimension>::ConstNodeListIterator nodeListItr = dataBase.nodeListBegin();
       nodeListItr != dataBase.nodeListEnd();
       ++nodeListItr) {
    controlIDs[*nodeListItr] = (*nodeListItr)->neighbor().coarseNeighborList();
    CHECK(controlIDs[*nodeListItr].size() == (*nodeListItr)->neighbor().numCoarse());
  }
  CHECK(controlIDs.size() == dataBase.numNodeLists());

  // Perform the check.
  return testNodeIteratorPair(dataBase.coarseNodeBegin(),
                              dataBase.coarseNodeEnd(),
                              controlIDs);
}

//------------------------------------------------------------------------------
// Test global RefineNodeIterators.
//------------------------------------------------------------------------------
template<typename Dimension>
std::string
testGlobalRefineNodeIterators(const DataBase<Dimension>& dataBase) {

  // Build the set of node IDs we expect the iterators to span.
  std::map<const NodeList<Dimension>*, std::vector<int> > controlIDs;
  for (typename DataBase<Dimension>::ConstNodeListIterator nodeListItr = dataBase.nodeListBegin();
       nodeListItr != dataBase.nodeListEnd();
       ++nodeListItr) {
    controlIDs[*nodeListItr] = (*nodeListItr)->neighbor().refineNeighborList();
    CHECK(controlIDs[*nodeListItr].size() == (*nodeListItr)->neighbor().numRefine());
  }
  CHECK(controlIDs.size() == dataBase.numNodeLists());

  // Perform the check.
  return testNodeIteratorPair(dataBase.refineNodeBegin(),
                              dataBase.refineNodeEnd(),
                              controlIDs);
}

}
