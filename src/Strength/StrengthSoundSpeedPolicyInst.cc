//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "StrengthSoundSpeedPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class StrengthSoundSpeedPolicy<Dim<1> >;
  template class StrengthSoundSpeedPolicy<Dim<2> >;
  template class StrengthSoundSpeedPolicy<Dim<3> >;
}

