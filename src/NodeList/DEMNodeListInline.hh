#include "Field/Field.hh"
//#include "Utilities/SpheralFunctions.hh"
//#include "SmoothingScaleBase.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Angular Velocity per node.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
Field<Dimension, typename Dimension::Vector>&
DEMNodeList<Dimension>::angularVelocity() {
  REQUIRE(mAngularVelocity.nodeListPtr() == this);
  return mAngularVelocity;
}

template<typename Dimension>
inline
const Field<Dimension, typename Dimension::Vector>&
DEMNodeList<Dimension>::angularVelocity() const {
  REQUIRE(mAngularVelocity.nodeListPtr() == this);
  return mAngularVelocity;
}

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

}
