//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Strength/YieldStrengthPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class YieldStrengthPolicy<Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class YieldStrengthPolicy<Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class YieldStrengthPolicy<Dim<3> >;
#endif
}