#include <algorithm>

namespace Spheral {

//------------------------------------------------------------------------------
// Increment the iterator.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
InternalNodeIterator<Dimension>&
InternalNodeIterator<Dimension>::
operator++() {

  // Increment the current node and cache indicies.
  ++mNodeID;

  // If the index is out of range, then proceed to the next NodeList with a 
  // nonzero number of nodes.
  if (mNodeID >= (int)(*mNodeListItr)->numInternalNodes()) {
    mNodeID = 0;
    ++mNodeListItr;
    while (mNodeListItr < mNodeListEnd &&
           (*mNodeListItr)->numInternalNodes() == 0) ++mNodeListItr;
    mFieldID = std::distance(mNodeListBegin, mNodeListItr);
  }
  
  ENSURE(valid());
  return *this;
}

}

