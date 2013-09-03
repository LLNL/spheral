//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "TensorDamagePolicy.cc"

namespace Spheral {
  template class TensorDamagePolicy<Dim<1> >;
  template class TensorDamagePolicy<Dim<2> >;
  template class TensorDamagePolicy<Dim<3> >;
}
