//---------------------------------Spheral++----------------------------------//
// RefineNodeIterator -- The version of the NodeIterator that goes over all
// refine nodes in a list of NodeLists based on their Neighbor states.
//
// Created by J. Michael Owen, Tue Mar 18 21:32:13 PST 2003
//----------------------------------------------------------------------------//
#include "RefineNodeIterator.hh"

#include <algorithm>
using std::vector;

namespace Spheral {

//------------------------------------------------------------------------------
// Default constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
RefineNodeIterator<Dimension>::
RefineNodeIterator():
  NodeIteratorBase<Dimension>(),
  mNodeIDItr(),
  mRefineNeighbors() {
}

//------------------------------------------------------------------------------
// Construct with the given NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
RefineNodeIterator<Dimension>::
RefineNodeIterator(typename vector<NodeList<Dimension>*>::const_iterator nodeListItr,
                   typename vector<NodeList<Dimension>*>::const_iterator nodeListBegin,
                   typename vector<NodeList<Dimension>*>::const_iterator nodeListEnd,
                   const std::vector<std::vector<int>>& refineNeighbors):
  NodeIteratorBase<Dimension>(),
  mNodeIDItr(),
  mRefineNeighbors(refineNeighbors) {
  initialize(nodeListItr,
             nodeListBegin,
             nodeListEnd,
             vector<int>::const_iterator(),
             refineNeighbors);
  ENSURE(valid());
}

//------------------------------------------------------------------------------
// Construct with the given NodeList and node ID.
//------------------------------------------------------------------------------
template<typename Dimension>
RefineNodeIterator<Dimension>::
RefineNodeIterator(typename vector<NodeList<Dimension>*>::const_iterator nodeListItr, 
                   typename vector<NodeList<Dimension>*>::const_iterator nodeListBegin,
                   typename vector<NodeList<Dimension>*>::const_iterator nodeListEnd,
                   vector<int>::const_iterator IDItr,
                   const std::vector<std::vector<int>>& refineNeighbors):
  NodeIteratorBase<Dimension>(),
  mNodeIDItr(),
  mRefineNeighbors(refineNeighbors) {
  initialize(nodeListItr,
             nodeListBegin,
             nodeListEnd,
             IDItr,
             refineNeighbors);
  ENSURE(valid());
}

//------------------------------------------------------------------------------
// Copy constructor
//------------------------------------------------------------------------------
template<typename Dimension>
RefineNodeIterator<Dimension>::
RefineNodeIterator(const RefineNodeIterator<Dimension>& rhs):
  NodeIteratorBase<Dimension>(rhs),
  mNodeIDItr(),
  mRefineNeighbors(rhs.mRefineNeighbors) {
  initialize(mNodeListItr,
             mNodeListBegin,
             mNodeListEnd,
             rhs.mNodeIDItr,
             rhs.mRefineNeighbors);
  ENSURE(valid() == rhs.valid());
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
RefineNodeIterator<Dimension>::
~RefineNodeIterator() {
}

//------------------------------------------------------------------------------
// Valid test.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
RefineNodeIterator<Dimension>::
valid() const {

  // NodeIteratorBase test.
  const bool baseTest = NodeIteratorBase<Dimension>::valid();

  // Verify that the node ID corresponds to a refine node.
  bool refineTest;
  if (mNodeListItr != mNodeListEnd) {
    refineTest = (mNodeID == *mNodeIDItr && 
                  find(mRefineNeighbors[mFieldID].begin(),
                       mRefineNeighbors[mFieldID].end(),
                       *mNodeIDItr) != mRefineNeighbors[mFieldID].end());
  } else {
    refineTest = mNodeID == 0;
  }

  return baseTest && refineTest;
}

//------------------------------------------------------------------------------
// Private initialization method.
//------------------------------------------------------------------------------
template<typename Dimension>
void
RefineNodeIterator<Dimension>::
initialize(typename vector<NodeList<Dimension>*>::const_iterator nodeListItr,
           typename vector<NodeList<Dimension>*>::const_iterator nodeListBegin,
           typename vector<NodeList<Dimension>*>::const_iterator nodeListEnd,
           vector<int>::const_iterator IDItr,
           const std::vector<std::vector<int>>& refineNeighbors) {

  // Pre-conditions.
  mFieldID = distance(nodeListBegin, nodeListItr);
  REQUIRE(nodeListItr == nodeListEnd ||
          (nodeListItr < nodeListEnd && 
           IDItr >= refineNeighbors[mFieldID].begin() &&
           IDItr <= refineNeighbors[mFieldID].end()));

  mNodeListBegin = nodeListBegin;
  mNodeListEnd = nodeListEnd;
  mNodeListItr = nodeListItr;
  if (nodeListItr < nodeListEnd) {
    mNodeIDItr = mRefineNeighbors[mFieldID].begin() + std::distance(refineNeighbors[mFieldID].begin(), IDItr);
    mNodeID = *IDItr;
  } else {
    mNodeIDItr = vector<int>::const_iterator();
    mNodeID = 0;
  }

  // Post conditions.
  ENSURE(valid());
}

}

