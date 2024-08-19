//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "RK/ContinuityVolumePolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class ContinuityVolumePolicy<Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class ContinuityVolumePolicy<Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class ContinuityVolumePolicy<Dim<3> >;
#endif
}