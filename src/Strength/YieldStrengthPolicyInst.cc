//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "YieldStrengthPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class YieldStrengthPolicy<Dim<1> >;
  template class YieldStrengthPolicy<Dim<2> >;
  template class YieldStrengthPolicy<Dim<3> >;
}

