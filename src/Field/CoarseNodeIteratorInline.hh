#include <algorithm>

namespace Spheral {

//------------------------------------------------------------------------------
// Increment the iterator
//------------------------------------------------------------------------------
template<typename Dimension>
inline
CoarseNodeIterator<Dimension>&
CoarseNodeIterator<Dimension>::operator++() {

  // Increment the current coarse node iterator.
  ++mNodeIDItr;

  // If the index is out of range, then proceed to the next NodeList with a 
  // nonzero number of internal nodes.
  if (mNodeListItr < mNodeListEnd && mNodeIDItr < mCoarseNeighbors[mFieldID].end()) {
    mNodeID = *mNodeIDItr;
  } else {
    ++mNodeListItr;
    ++mFieldID;
    while (mNodeListItr < mNodeListEnd &&
           mCoarseNeighbors[mFieldID].size() == 0) {
      ++mNodeListItr;
      ++mFieldID;
    }
    if (mNodeListItr < mNodeListEnd) {
      mNodeIDItr = mCoarseNeighbors[mFieldID].begin();
      if (mNodeIDItr < mCoarseNeighbors[mFieldID].end()) {
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
