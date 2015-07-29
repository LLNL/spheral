//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "HullVolumePolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class HullVolumePolicy<Dim<1> >;
  template class HullVolumePolicy<Dim<2> >;
  template class HullVolumePolicy<Dim<3> >;
}

