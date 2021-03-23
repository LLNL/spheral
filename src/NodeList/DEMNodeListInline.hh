#include "Field/Field.hh"
#include "Utilities/SpheralFunctions.hh"
#include "SmoothingScaleBase.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Mass density per node.
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

}
