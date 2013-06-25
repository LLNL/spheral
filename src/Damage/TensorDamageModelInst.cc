//------------------------------------------------------------------------------
// Instantiations.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "TensorDamageModel.cc"

namespace Spheral {
  template class PhysicsSpace::TensorDamageModel<Dim<1> >;
  template class PhysicsSpace::TensorDamageModel<Dim<2> >;
  template class PhysicsSpace::TensorDamageModel<Dim<3> >;
}
