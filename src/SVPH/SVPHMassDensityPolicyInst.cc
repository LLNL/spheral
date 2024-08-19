//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "SVPH/SVPHMassDensityPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class SVPHMassDensityPolicy<Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class SVPHMassDensityPolicy<Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class SVPHMassDensityPolicy<Dim<3> >;
#endif
}