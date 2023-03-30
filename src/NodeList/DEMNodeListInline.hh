#include "Field/Field.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// particle radius per node.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
Field<Dimension, typename Dimension::Scalar>&
DEMNodeList<Dimension>::particleRadius() {
  REQUIRE(mParticleRadius.nodeListPtr() == this);
  return mParticleRadius;
}

template<typename Dimension>
inline
const Field<Dimension, typename Dimension::Scalar>&
DEMNodeList<Dimension>::particleRadius() const {
  REQUIRE(mParticleRadius.nodeListPtr() == this);
  return mParticleRadius;
}


//------------------------------------------------------------------------------
// composite particle indices
//------------------------------------------------------------------------------
template<typename Dimension>
inline
Field<Dimension, int>&
DEMNodeList<Dimension>::compositeParticleIndex() {
  REQUIRE(mCompositeParticleIndex.nodeListPtr() == this);
  return mCompositeParticleIndex;
}

template<typename Dimension>
inline
const Field<Dimension, int>&
DEMNodeList<Dimension>::compositeParticleIndex() const {
  REQUIRE(mCompositeParticleIndex.nodeListPtr() == this);
  return mCompositeParticleIndex;
}
}
