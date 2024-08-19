//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Damage/YoungsModulusPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class YoungsModulusPolicy<Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class YoungsModulusPolicy<Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class YoungsModulusPolicy<Dim<3> >;
#endif
}