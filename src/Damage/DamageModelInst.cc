//------------------------------------------------------------------------------
// Instantiations.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "DamageModel.cc"

namespace Spheral {
  template class PhysicsSpace::DamageModel<Dim<1> >;
  template class PhysicsSpace::DamageModel<Dim<2> >;
  template class PhysicsSpace::DamageModel<Dim<3> >;
}
