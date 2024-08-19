//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "RK/HVolumePolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class HVolumePolicy<Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class HVolumePolicy<Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class HVolumePolicy<Dim<3> >;
#endif
}