//---------------------------------Spheral++----------------------------------//
// CoarseNodeIterator -- The version of the NodeIterator that goes over all
// coarse nodes in a list of NodeLists based on their Neighbor states.
//
// Created by J. Michael Owen, Tue Mar 18 21:32:13 PST 2003
//----------------------------------------------------------------------------//
#include <algorithm>

#include "CoarseNodeIterator.hh"

namespace Spheral {

using namespace std;
using NodeSpace::NodeList;

//------------------------------------------------------------------------------
// Default constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
CoarseNodeIterator<Dimension>::
CoarseNodeIterator():
  NodeIteratorBase<Dimension>(),
  mNodeIDItr() {
}

//------------------------------------------------------------------------------
// Construct with the given NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
CoarseNodeIterator<Dimension>::
CoarseNodeIterator(typename vector<NodeList<Dimension>*>::const_iterator nodeListItr,
                   typename vector<NodeList<Dimension>*>::const_iterator nodeListBegin,
                   typename vector<NodeList<Dimension>*>::const_iterator nodeListEnd):
  NodeIteratorBase<Dimension>(),
  mNodeIDItr() {
  initialize(nodeListItr,
             nodeListBegin,
             nodeListEnd,
             vector<int>::const_iterator());
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
                   vector<int>::const_iterator IDItr):
  NodeIteratorBase<Dimension>(),
  mNodeIDItr(IDItr) {
  initialize(nodeListItr,
             nodeListBegin,
             nodeListEnd,
             IDItr);
  ENSURE(valid());
}

//------------------------------------------------------------------------------
// Copy constructor
//------------------------------------------------------------------------------
template<typename Dimension>
CoarseNodeIterator<Dimension>::
CoarseNodeIterator(const CoarseNodeIterator<Dimension>& rhs):
  NodeIteratorBase<Dimension>(rhs),
  mNodeIDItr(rhs.mNodeIDItr) {
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
                  find((*mNodeListItr)->neighbor().coarseNeighborBegin(),
                       (*mNodeListItr)->neighbor().coarseNeighborEnd(),
                       *mNodeIDItr) != (*mNodeListItr)->neighbor().coarseNeighborEnd());
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
           vector<int>::const_iterator IDItr) {

  // Pre-conditions.
  REQUIRE((nodeListItr == nodeListEnd && IDItr == vector<int>::const_iterator()) ||
          (nodeListItr < nodeListEnd && 
           IDItr >= (*nodeListItr)->neighbor().coarseNeighborBegin() &&
           IDItr <= (*nodeListItr)->neighbor().coarseNeighborEnd()));

  if (nodeListItr < nodeListEnd) {
    mNodeID = *IDItr;
  } else {
    mNodeID = 0;
  }    
  mFieldID = distance(nodeListBegin, nodeListItr);
  mNodeListBegin = nodeListBegin;
  mNodeListEnd = nodeListEnd;
  mNodeListItr = nodeListItr;

  // Post conditions.
  ENSURE(valid());
}

}

