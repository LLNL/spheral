#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "NodeList/NodeList.hh"

namespace Spheral {
namespace PhysicsSpace {

//------------------------------------------------------------------------------
// Return the acceleration.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Vector
ConstantAcceleration<Dimension>::
a0() const {
  return ma0;
}

//------------------------------------------------------------------------------
// The NodeList this acceleration applies to.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const NodeSpace::NodeList<Dimension>&
ConstantAcceleration<Dimension>::
nodeList() const {
  REQUIRE(mNodeListPtr != 0);
  return *mNodeListPtr;
}

//------------------------------------------------------------------------------
// Flags for the set of nodes the acceleration is to be applied to.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldSpace::Field<Dimension, int>&
ConstantAcceleration<Dimension>::
flags() const {
  return mFlags;
}

}
}
