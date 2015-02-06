//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "HVolumePolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class HVolumePolicy<Dim<1> >;
  template class HVolumePolicy<Dim<2> >;
  template class HVolumePolicy<Dim<3> >;
}

