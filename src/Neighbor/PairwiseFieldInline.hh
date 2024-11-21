//---------------------------------Spheral++----------------------------------//
// PairwiseField
//
// Stores a value per node pair in a NodePairList. Because connectivity is
// allowed to change step to step in our meshfree methods, PairwiseField is
// ephemeral and will be invalidated when topology is updated.
//
// Created by J. Michael Owen, Wed Nov 20 14:44:44 PST 2024
//----------------------------------------------------------------------------//

#include "Neighbor/ConnectivityMap.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
PairwiseField<Dimension, Value>::PairwiseField(const ConnectivityMap<Dimension>& connectivity):
  mPairsPtr(connectivity.nodePairListPtr()),
  mValues() {
  if (auto p = mPairsPtr.lock()) {
    mValues.resize(p->size());
  } else {
    VERIFY2(false, "PairwiseField constructed with invalid NodePairList");
  }
}

//------------------------------------------------------------------------------
// Index by pair
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
inline
const Value&
PairwiseField<Dimension, Value>::operator()(const NodePairIdxType& x) const {
  if (auto p = mPairsPtr.lock()) {
    return mValues[p->index(x)];
  }
  VERIFY2(false, "PairwiseField ERROR: attempt to index with invalid pair " << x);
}

template<typename Dimension, typename Value>
inline
Value&
PairwiseField<Dimension, Value>::operator()(const NodePairIdxType& x) {
  if (auto p = mPairsPtr.lock()) {
    return mValues[p->index(x)];
  }
  VERIFY2(false, "PairwiseField ERROR: attempt to index with invalid pair " << x);
}

//------------------------------------------------------------------------------
// NodePairList
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
inline
const NodePairList&
PairwiseField<Dimension, Value>::pairs() const {
  if (auto p = mPairsPtr.lock()) {
    return *p;
  }
  VERIFY2(false, "Orphaned PairwiseField without NodePairList");
}

}
