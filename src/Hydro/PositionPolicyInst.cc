//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "PositionPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class PositionPolicy<Dim<1> >;
  template class PositionPolicy<Dim<2> >;
  template class PositionPolicy<Dim<3> >;
}
