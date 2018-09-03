//---------------------------------Spheral++----------------------------------//
// CoarseNodeIterator -- The version of the NodeIterator that goes over all
// coarse nodes in a list of NodeLists based on their Neighbor states.
//
// Created by J. Michael Owen, Tue Mar 18 21:32:13 PST 2003
//----------------------------------------------------------------------------//
#include "CoarseNodeIterator.hh"

#include <algorithm>
using std::vector;

namespace Spheral {

//------------------------------------------------------------------------------
// Default constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
CoarseNodeIterator<Dimension>::
CoarseNodeIterator():
  NodeIteratorBase<Dimension>(),
  mNodeIDItr(),
  mCoarseNeighbors() {
}

//------------------------------------------------------------------------------
// Construct with the given NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
CoarseNodeIterator<Dimension>::
CoarseNodeIterator(typename vector<NodeList<Dimension>*>::const_iterator nodeListItr,
                   typename vector<NodeList<Dimension>*>::const_iterator nodeListBegin,
                   typename vector<NodeList<Dimension>*>::const_iterator nodeListEnd,
                   const std::vector<std::vector<int>>& coarseNeighbors):
  NodeIteratorBase<Dimension>(),
  mNodeIDItr(),
  mCoarseNeighbors(coarseNeighbors) {
  initialize(nodeListItr,
             nodeListBegin,
             nodeListEnd,
             vector<int>::const_iterator(),
             coarseNeighbors);
  ENSURE(valid());
}

//------------------------------------------------------------------------------
// Construct with the given NodeList and node ID.
//------------------------------------------------------------------------------
template<typename Dimension>
CoarseNodeIterator<Dimension>::
CoarseNodeIterator(typename vector<NodeList<Dimension>*>::const_iterator nodeListItr, 
                   typename vector<NodeList<Dimension>*>::const_iterator nodeListBegin,
                   typename vector<NodeList<Dimension>*>::const_iterator nodeListEnd,
                   vector<int>::const_iterator IDItr,
                   const std::vector<std::vector<int>>& coarseNeighbors):
  NodeIteratorBase<Dimension>(),
  mNodeIDItr(),
  mCoarseNeighbors(coarseNeighbors) {
  initialize(nodeListItr,
             nodeListBegin,
             nodeListEnd,
             IDItr,
             coarseNeighbors);
  ENSURE(valid());
}

//------------------------------------------------------------------------------
// Copy constructor
//------------------------------------------------------------------------------
template<typename Dimension>
CoarseNodeIterator<Dimension>::
CoarseNodeIterator(const CoarseNodeIterator<Dimension>& rhs):
  NodeIteratorBase<Dimension>(rhs),
  mNodeIDItr(),
  mCoarseNeighbors(rhs.mCoarseNeighbors) {
  initialize(mNodeListItr,
             mNodeListBegin,
             mNodeListEnd,
             rhs.mNodeIDItr,
             rhs.mCoarseNeighbors);
  ENSURE(valid() == rhs.valid());
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
CoarseNodeIterator<Dimension>::
~CoarseNodeIterator() {
}

//------------------------------------------------------------------------------
// Valid test.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
CoarseNodeIterator<Dimension>::
valid() const {

  // NodeIteratorBase test.
  const bool baseTest = NodeIteratorBase<Dimension>::valid();

  // Verify that the node ID corresponds to a coarse node.
  bool coarseTest;
  if (mNodeListItr != mNodeListEnd) {
    coarseTest = (mNodeID == *mNodeIDItr && 
                  find(mCoarseNeighbors[mFieldID].begin(),
                       mCoarseNeighbors[mFieldID].end(),
                       *mNodeIDItr) != mCoarseNeighbors[mFieldID].end());
  } else {
    coarseTest = mNodeID == 0;
  }

  return baseTest && coarseTest;
}

//------------------------------------------------------------------------------
// Private initialization method.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CoarseNodeIterator<Dimension>::
initialize(typename vector<NodeList<Dimension>*>::const_iterator nodeListItr,
           typename vector<NodeList<Dimension>*>::const_iterator nodeListBegin,
           typename vector<NodeList<Dimension>*>::const_iterator nodeListEnd,
           vector<int>::const_iterator IDItr,
           const std::vector<std::vector<int>>& coarseNeighbors) {

  // Pre-conditions.
  mFieldID = distance(nodeListBegin, nodeListItr);
  REQUIRE(nodeListItr == nodeListEnd ||
          (nodeListItr < nodeListEnd && 
           IDItr >= coarseNeighbors[mFieldID].begin() &&
           IDItr <= coarseNeighbors[mFieldID].end()));

  mNodeListBegin = nodeListBegin;
  mNodeListEnd = nodeListEnd;
  mNodeListItr = nodeListItr;
  if (nodeListItr < nodeListEnd) {
    mNodeIDItr = mCoarseNeighbors[mFieldID].begin() + std::distance(coarseNeighbors[mFieldID].begin(), IDItr);
    mNodeID = *IDItr;
  } else {
    mNodeIDItr = vector<int>::const_iterator();
    mNodeID = 0;
  }

  // Post conditions.
  ENSURE(valid());
}

}

