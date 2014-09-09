//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "PressurePolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class PressurePolicy<Dim<1> >;
  template class PressurePolicy<Dim<2> >;
  template class PressurePolicy<Dim<3> >;
}
