#include <algorithm>

namespace Spheral {

//------------------------------------------------------------------------------
// Increment the iterator
//------------------------------------------------------------------------------
template<typename Dimension>
inline
RefineNodeIterator<Dimension>&
RefineNodeIterator<Dimension>::operator++() {

  // Increment the current refine node iterator.
  ++mNodeIDItr;

  // If the index is out of range, then proceed to the next NodeList with a 
  // nonzero number of internal nodes.
  if (mNodeListItr < mNodeListEnd && mNodeIDItr < mRefineNeighbors[mFieldID].end()) {
    mNodeID = *mNodeIDItr;
  } else {
    ++mNodeListItr;
    ++mFieldID;
    while (mNodeListItr < mNodeListEnd &&
           mRefineNeighbors[mFieldID].size() == 0) {
      ++mNodeListItr;
      ++mFieldID;
    }
    if (mNodeListItr < mNodeListEnd) {
      mNodeIDItr = mRefineNeighbors[mFieldID].begin();
      if (mNodeIDItr < mRefineNeighbors[mFieldID].end()) {
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

