//---------------------------------Spheral++----------------------------------//
// InternalNodeIterator -- The version of the NodeIterator that goes over all nodes
// in a list of NodeLists.
//
// Created by J. Michael Owen, Mon Mar 17 23:13:08 PST 2003
//----------------------------------------------------------------------------//
#include "InternalNodeIterator.hh"

#include <algorithm>
using std::vector;

namespace Spheral {

//------------------------------------------------------------------------------
// Default constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
InternalNodeIterator<Dimension>::InternalNodeIterator():
  NodeIteratorBase<Dimension>() {
}

//------------------------------------------------------------------------------
// Construct with the given NodeList and node ID.
//------------------------------------------------------------------------------
template<typename Dimension>
InternalNodeIterator<Dimension>::
InternalNodeIterator(typename vector<NodeList<Dimension>*>::const_iterator nodeListItr, 
                     typename vector<NodeList<Dimension>*>::const_iterator nodeListBegin,
                     typename vector<NodeList<Dimension>*>::const_iterator nodeListEnd,
                     int nodeID):
  NodeIteratorBase<Dimension>() {

  // Pre-conditions.
  REQUIRE((nodeListItr == nodeListEnd && nodeID == 0) ||
          (nodeListItr < nodeListEnd &&
           nodeID >= 0 && nodeID <= (int)(*nodeListItr)->numInternalNodes()));

  mNodeID = nodeID;
  mFieldID = distance(nodeListBegin, nodeListItr);
  mNodeListBegin = nodeListBegin;
  mNodeListEnd = nodeListEnd;
  mNodeListItr = nodeListItr;

  // Post conditions.
  ENSURE(valid());
}

//------------------------------------------------------------------------------
// Copy constructor
//------------------------------------------------------------------------------
template<typename Dimension>
InternalNodeIterator<Dimension>::
InternalNodeIterator(const InternalNodeIterator<Dimension>& rhs):
  NodeIteratorBase<Dimension>(rhs) {
  ENSURE(valid() == rhs.valid());
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
InternalNodeIterator<Dimension>::
~InternalNodeIterator() {
}

//------------------------------------------------------------------------------
// Valid test.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
InternalNodeIterator<Dimension>::
valid() const {

  // NodeIteratorBase test.
  const bool baseTest = NodeIteratorBase<Dimension>::valid();

  // Verify that the node ID corresponds to an internal node.
  bool internalTest;
  if (mNodeListItr != mNodeListEnd) {
    internalTest = this->internalNode();
  } else {
    internalTest = mNodeID == 0;
  }

  return baseTest && internalTest;
}

}

