#include <algorithm>

namespace Spheral {

//------------------------------------------------------------------------------
// Increment the iterator.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
GhostNodeIterator<Dimension>&
GhostNodeIterator<Dimension>::
operator++() {

  // Increment the current node and cache indicies.
  ++mNodeID;

  // If the index is out of range, then proceed to the next NodeList with a 
  // nonzero number of nodes.
  if (mNodeID >= (int)(*mNodeListItr)->numNodes()) {
    ++mNodeListItr;
    while (mNodeListItr < mNodeListEnd &&
           (*mNodeListItr)->numGhostNodes() == 0) ++mNodeListItr;
    mFieldID = std::distance(mNodeListBegin, mNodeListItr);
    if (mNodeListItr < mNodeListEnd) {
      mNodeID = (*mNodeListItr)->firstGhostNode();
    } else {
      mNodeID = 0;
    }
  }
  
  ENSURE(valid());
  return *this;
}

}
