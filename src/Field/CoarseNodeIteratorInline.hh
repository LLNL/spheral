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
  if (mNodeListItr < mNodeListEnd && mNodeIDItr < (*mNodeListItr)->neighbor().coarseNeighborEnd()) {
    mNodeID = *mNodeIDItr;
  } else {
    ++mNodeListItr;
    while (mNodeListItr < mNodeListEnd &&
           (*mNodeListItr)->neighbor().numCoarse() == 0) ++mNodeListItr;
    mFieldID = std::distance(mNodeListBegin, mNodeListItr);
    if (mNodeListItr < mNodeListEnd) {
      mNodeIDItr = (*mNodeListItr)->neighbor().coarseNeighborBegin();
      if (mNodeIDItr < (*mNodeListItr)->neighbor().coarseNeighborEnd()) {
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
