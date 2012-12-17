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
  if (mNodeListItr < mNodeListEnd && mNodeIDItr < (*mNodeListItr)->neighbor().refineNeighborEnd()) {
    mNodeID = *mNodeIDItr;
  } else {
    ++mNodeListItr;
    while (mNodeListItr < mNodeListEnd &&
           (*mNodeListItr)->neighbor().numRefine() == 0) ++mNodeListItr;
    mFieldID = std::distance(mNodeListBegin, mNodeListItr);
    if (mNodeListItr < mNodeListEnd) {
      mNodeIDItr = (*mNodeListItr)->neighbor().refineNeighborBegin();
      if (mNodeIDItr < (*mNodeListItr)->neighbor().refineNeighborEnd()) {
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

