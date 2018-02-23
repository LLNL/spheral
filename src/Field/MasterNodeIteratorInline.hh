#include <algorithm>

namespace Spheral {

//------------------------------------------------------------------------------
// Increment the iterator
//------------------------------------------------------------------------------
template<typename Dimension>
inline
MasterNodeIterator<Dimension>&
MasterNodeIterator<Dimension>::operator++() {

  // Increment the current master node iterator.
  ++mNodeIDItr;

  // If the index is out of range, then proceed to the next NodeList with a 
  // nonzero number of internal nodes.
  if (mNodeListItr < mNodeListEnd && mNodeIDItr < mMasterLists[mFieldID].end()) {
    mNodeID = *mNodeIDItr;
  } else {
    ++mNodeListItr;
    ++mFieldID;
    while (mNodeListItr < mNodeListEnd &&
           mMasterLists[mFieldID].size() == 0) {
      ++mNodeListItr;
      ++mFieldID;
    }
    if (mNodeListItr < mNodeListEnd) {
      mNodeIDItr = mMasterLists[mFieldID].begin();
      if (mNodeIDItr < mMasterLists[mFieldID].end()) {
        mNodeID = *mNodeIDItr;
      } else {
        mNodeID = 0;
      }
    } else {
      mNodeListItr = mNodeListEnd;
      mNodeID = 0;
    }
  }
  
  ENSURE(valid());
  return *this;
}

}

