//---------------------------------Spheral++----------------------------------//
// MasterNodeIterator -- The version of the NodeIterator that goes over all
// master nodes in a list of NodeLists based on their Neighbor states.
//
// Created by J. Michael Owen, Mon Mar 17 23:13:08 PST 2003
//----------------------------------------------------------------------------//
#include "MasterNodeIterator.hh"

#include <algorithm>
using std::vector;

namespace Spheral {

//------------------------------------------------------------------------------
// Default constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
MasterNodeIterator<Dimension>::
MasterNodeIterator():
  NodeIteratorBase<Dimension>(),
  mNodeIDItr(),
  mMasterLists() {
}

//------------------------------------------------------------------------------
// Construct with the given NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
MasterNodeIterator<Dimension>::
MasterNodeIterator(typename vector<NodeList<Dimension>*>::const_iterator nodeListItr,
                   typename vector<NodeList<Dimension>*>::const_iterator nodeListBegin,
                   typename vector<NodeList<Dimension>*>::const_iterator nodeListEnd,
                   const std::vector<std::vector<int>>& masterLists):
  NodeIteratorBase<Dimension>(),
  mNodeIDItr(),
  mMasterLists(masterLists) {
  initialize(nodeListItr,
             nodeListBegin,
             nodeListEnd,
             vector<int>::const_iterator(),
             masterLists);
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
                   vector<int>::const_iterator IDItr,
                   const std::vector<std::vector<int>>& masterLists):
  NodeIteratorBase<Dimension>(),
  mNodeIDItr(),
  mMasterLists(masterLists) {
  initialize(nodeListItr,
             nodeListBegin,
             nodeListEnd,
             IDItr,
             masterLists);
  ENSURE(valid());
}

//------------------------------------------------------------------------------
// Copy constructor
//------------------------------------------------------------------------------
template<typename Dimension>
MasterNodeIterator<Dimension>::
MasterNodeIterator(const MasterNodeIterator<Dimension>& rhs):
  NodeIteratorBase<Dimension>(rhs),
  mNodeIDItr(),
  mMasterLists(rhs.mMasterLists) {
  initialize(mNodeListItr,
             mNodeListBegin,
             mNodeListEnd,
             rhs.mNodeIDItr,
             rhs.mMasterLists);
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
                  find(mMasterLists[mFieldID].begin(),
                       mMasterLists[mFieldID].end(),
                       *mNodeIDItr) != mMasterLists[mFieldID].end());
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
           vector<int>::const_iterator IDItr,
           const std::vector<std::vector<int>>& masterLists) {

  // Pre-conditions.
  mFieldID = distance(nodeListBegin, nodeListItr);
  REQUIRE(nodeListItr == nodeListEnd ||
          (nodeListItr < nodeListEnd && 
           IDItr >= masterLists[mFieldID].begin() &&
           IDItr <= masterLists[mFieldID].end()));

  mNodeListBegin = nodeListBegin;
  mNodeListEnd = nodeListEnd;
  mNodeListItr = nodeListItr;
  if (nodeListItr < nodeListEnd) {
    mNodeIDItr = mMasterLists[mFieldID].begin() + std::distance(masterLists[mFieldID].begin(), IDItr);
    mNodeID = *IDItr;
  } else {
    mNodeIDItr = vector<int>::const_iterator();
    mNodeID = 0;
  }

  // Post conditions.
  ENSURE(valid());
}

}

