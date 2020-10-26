//---------------------------------Spheral++----------------------------------//
// AllNodeIterator -- The version of the NodeIterator that goes over all nodes
// in a list of NodeLists.
//
// Created by J. Michael Owen, Mon Mar 17 17:31:14 PST 2003
//----------------------------------------------------------------------------//
#include "AllNodeIterator.hh"

#include <algorithm>
using std::vector;

namespace Spheral {

//------------------------------------------------------------------------------
// Default constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
AllNodeIterator<Dimension>::AllNodeIterator():
  NodeIteratorBase<Dimension>() {
}

//------------------------------------------------------------------------------
// Construct with the given NodeList and node ID.
//------------------------------------------------------------------------------
template<typename Dimension>
AllNodeIterator<Dimension>::
AllNodeIterator(typename vector<NodeList<Dimension>*>::const_iterator nodeListItr, 
                typename vector<NodeList<Dimension>*>::const_iterator nodeListBegin,
                typename vector<NodeList<Dimension>*>::const_iterator nodeListEnd,
                int nodeID):
  NodeIteratorBase<Dimension>() {

  // Pre-conditions.
  REQUIRE((nodeListItr == nodeListEnd && nodeID == 0) ||
          (nodeListItr < nodeListEnd && 
           nodeID >= 0 && nodeID <= (int)(*nodeListItr)->numNodes()));

  mNodeID = nodeID;
  mFieldID = distance(nodeListBegin, nodeListItr);
  mNodeListBegin = nodeListBegin;
  mNodeListEnd = nodeListEnd;
  mNodeListItr = nodeListItr;
  
  // Post conditions.
  ENSURE(this->valid());
}

//------------------------------------------------------------------------------
// Copy constructor
//------------------------------------------------------------------------------
template<typename Dimension>
AllNodeIterator<Dimension>::
AllNodeIterator(const AllNodeIterator<Dimension>& rhs):
  NodeIteratorBase<Dimension>(rhs) {
  ENSURE(this->valid() == rhs.valid());
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
AllNodeIterator<Dimension>::
~AllNodeIterator() {
}

}

