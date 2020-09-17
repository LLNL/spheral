#include <algorithm>

namespace Spheral {

//------------------------------------------------------------------------------
// Increment the iterator.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
AllNodeIterator<Dimension>&
AllNodeIterator<Dimension>::
operator++() {

  // Increment the current node and cache indicies.
  ++mNodeID;

  // If the index is out of range, then proceed to the next NodeList with a 
  // nonzero number of nodes.
  if (mNodeID >= (int)(*mNodeListItr)->numNodes()) {
    mNodeID = 0;
    ++mNodeListItr;
    while (mNodeListItr < mNodeListEnd &&
           (*mNodeListItr)->numNodes() == 0) ++mNodeListItr;
    mFieldID = std::distance(mNodeListBegin, mNodeListItr);
  }
  
  ENSURE(this->valid());
  return *this;
}

}

