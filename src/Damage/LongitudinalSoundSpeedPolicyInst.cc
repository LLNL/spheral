//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "LongitudinalSoundSpeedPolicy.cc"

namespace Spheral {
  template class LongitudinalSoundSpeedPolicy<Dim<1> >;
  template class LongitudinalSoundSpeedPolicy<Dim<2> >;
  template class LongitudinalSoundSpeedPolicy<Dim<3> >;
}
