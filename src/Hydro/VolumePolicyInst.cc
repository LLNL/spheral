//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "VolumePolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class VolumePolicy<Dim<1> >;
  template class VolumePolicy<Dim<2> >;
  template class VolumePolicy<Dim<3> >;
}

