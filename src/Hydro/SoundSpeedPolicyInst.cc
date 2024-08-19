//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Hydro/SoundSpeedPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class SoundSpeedPolicy<Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class SoundSpeedPolicy<Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class SoundSpeedPolicy<Dim<3> >;
#endif
}