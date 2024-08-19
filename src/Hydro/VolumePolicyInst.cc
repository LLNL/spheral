//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Hydro/VolumePolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class VolumePolicy<Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class VolumePolicy<Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class VolumePolicy<Dim<3> >;
#endif
}