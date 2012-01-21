//---------------------------------Spheral++----------------------------------//
// GhostNodeIterator -- The version of the NodeIterator that goes over all nodes
// in a list of NodeLists.
//
// Created by J. Michael Owen, Mon Mar 17 17:31:14 PST 2003
//----------------------------------------------------------------------------//
#include <algorithm>

#include "GhostNodeIterator.hh"

namespace Spheral {

using namespace std;
using NodeSpace::NodeList;

//------------------------------------------------------------------------------
// Default constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
GhostNodeIterator<Dimension>::GhostNodeIterator():
  NodeIteratorBase<Dimension>() {
}

//------------------------------------------------------------------------------
// Construct with the given NodeList and node ID.
//------------------------------------------------------------------------------
template<typename Dimension>
GhostNodeIterator<Dimension>::
GhostNodeIterator(typename vector<NodeList<Dimension>*>::const_iterator nodeListItr, 
                  typename vector<NodeList<Dimension>*>::const_iterator nodeListBegin,
                  typename vector<NodeList<Dimension>*>::const_iterator nodeListEnd,
                  int nodeID):
  NodeIteratorBase<Dimension>() {

  // Pre-conditions.
  REQUIRE((nodeListItr == nodeListEnd && nodeID == 0) ||
          (nodeListItr < nodeListEnd &&
           nodeID >= (*nodeListItr)->firstGhostNode() && nodeID < (*nodeListItr)->numNodes()));

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
GhostNodeIterator<Dimension>::
GhostNodeIterator(const GhostNodeIterator<Dimension>& rhs):
  NodeIteratorBase<Dimension>(rhs) {
  ENSURE(valid() == rhs.valid());
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
GhostNodeIterator<Dimension>::
~GhostNodeIterator() {
}

//------------------------------------------------------------------------------
// Valid test.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
GhostNodeIterator<Dimension>::
valid() const {

  // NodeIteratorBase test.
  const bool baseTest = NodeIteratorBase<Dimension>::valid();

  // Verify that the node ID corresponds to an ghost node.
  bool ghostTest;
  if (mNodeListItr != mNodeListEnd) {
    ghostTest = this->ghostNode();
  } else {
    ghostTest = mNodeID == 0;
  }

  return baseTest && ghostTest;
}

}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
namespace Spheral {
  template class GhostNodeIterator< Dim<1> >;
  template class GhostNodeIterator< Dim<2> >;
  template class GhostNodeIterator< Dim<3> >;
}
