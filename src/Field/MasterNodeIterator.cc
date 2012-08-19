//---------------------------------Spheral++----------------------------------//
// MasterNodeIterator -- The version of the NodeIterator that goes over all
// master nodes in a list of NodeLists based on their Neighbor states.
//
// Created by J. Michael Owen, Mon Mar 17 23:13:08 PST 2003
//----------------------------------------------------------------------------//
#include <algorithm>

#include "MasterNodeIterator.hh"

namespace Spheral {

using namespace std;
using NodeSpace::NodeList;

//------------------------------------------------------------------------------
// Default constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
MasterNodeIterator<Dimension>::
MasterNodeIterator():
  NodeIteratorBase<Dimension>(),
  mNodeIDItr(0) {
}

//------------------------------------------------------------------------------
// Construct with the given NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
MasterNodeIterator<Dimension>::
MasterNodeIterator(typename vector<NodeList<Dimension>*>::const_iterator nodeListItr,
                   typename vector<NodeList<Dimension>*>::const_iterator nodeListBegin,
                   typename vector<NodeList<Dimension>*>::const_iterator nodeListEnd):
  NodeIteratorBase<Dimension>(),
  mNodeIDItr(0) {
  initialize(nodeListItr,
             nodeListBegin,
             nodeListEnd,
             vector<int>::const_iterator(0));
  ENSURE(valid());
}

//------------------------------------------------------------------------------
// Construct with the given NodeList and node ID.
//------------------------------------------------------------------------------
template<typename Dimension>
MasterNodeIterator<Dimension>::
MasterNodeIterator(typename vector<NodeList<Dimension>*>::const_iterator nodeListItr, 
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
MasterNodeIterator<Dimension>::
MasterNodeIterator(const MasterNodeIterator<Dimension>& rhs):
  NodeIteratorBase<Dimension>(rhs),
  mNodeIDItr(rhs.mNodeIDItr) {
  ENSURE(valid() == rhs.valid());
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
MasterNodeIterator<Dimension>::
~MasterNodeIterator() {
}

//------------------------------------------------------------------------------
// Valid test.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
MasterNodeIterator<Dimension>::
valid() const {

  // NodeIteratorBase test.
  const bool baseTest = NodeIteratorBase<Dimension>::valid();

  // Verify that the node ID corresponds to a master node.
  bool masterTest;
  if (mNodeListItr != mNodeListEnd) {
    masterTest = (mNodeID == *mNodeIDItr && 
                  find((*mNodeListItr)->neighbor().masterBegin(),
                       (*mNodeListItr)->neighbor().masterEnd(),
                       *mNodeIDItr) != (*mNodeListItr)->neighbor().masterEnd());
  } else {
    masterTest = mNodeID == 0;
  }

  return baseTest && masterTest;
}

//------------------------------------------------------------------------------
// Private initialization method.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MasterNodeIterator<Dimension>::
initialize(typename vector<NodeList<Dimension>*>::const_iterator nodeListItr,
           typename vector<NodeList<Dimension>*>::const_iterator nodeListBegin,
           typename vector<NodeList<Dimension>*>::const_iterator nodeListEnd,
           vector<int>::const_iterator IDItr) {

  // Pre-conditions.
  REQUIRE((nodeListItr == nodeListEnd && IDItr == vector<int>::const_iterator(0)) ||
          (nodeListItr < nodeListEnd && 
           IDItr >= (*nodeListItr)->neighbor().masterBegin() &&
           IDItr <= (*nodeListItr)->neighbor().masterEnd()));

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

