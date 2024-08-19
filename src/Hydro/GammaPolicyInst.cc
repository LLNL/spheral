//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Hydro/GammaPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class GammaPolicy<Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class GammaPolicy<Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class GammaPolicy<Dim<3> >;
#endif
}