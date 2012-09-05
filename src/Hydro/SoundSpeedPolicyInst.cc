//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SoundSpeedPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class SoundSpeedPolicy<Dim<1> >;
  template class SoundSpeedPolicy<Dim<2> >;
  template class SoundSpeedPolicy<Dim<3> >;
}

