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
  if (mNodeListItr < mNodeListEnd && mNodeIDItr < (*mNodeListItr)->neighbor().masterEnd()) {
    mNodeID = *mNodeIDItr;
  } else {
    ++mNodeListItr;
    while (mNodeListItr < mNodeListEnd &&
           (*mNodeListItr)->neighbor().numMaster() == 0) ++mNodeListItr;
    mFieldID = std::distance(mNodeListBegin, mNodeListItr);
    if (mNodeListItr < mNodeListEnd) {
      mNodeIDItr = (*mNodeListItr)->neighbor().masterBegin();
      if (mNodeIDItr < (*mNodeListItr)->neighbor().masterEnd()) {
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

