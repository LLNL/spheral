#include "NodePairList.hh"
#include "Utilities/DBC.hh"

namespace Spheral {
  
//------------------------------------------------------------------------------
// index
//------------------------------------------------------------------------------
size_t
NodePairList::index(const NodePairIdxType& x) const {
  if (mPair2Index.size() != mNodePairList.size()) computeLookup();  // Lazy evaluation
  auto itr = mPair2Index.find(x);
  CHECK(itr != mPair2Index.end());
  return itr->second;
}

//------------------------------------------------------------------------------
// Recompute the lookup table for NodePair->index
//------------------------------------------------------------------------------
void
NodePairList::computeLookup() const {
  mPair2Index.clear();
  const auto n = this->size();
  for (size_t k = 0u; k < n; ++k) {
    mPair2Index[mNodePairList[k]] = k;
  }
}

}
