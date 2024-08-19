//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Damage/ProbabilisticDamagePolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class ProbabilisticDamagePolicy<Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class ProbabilisticDamagePolicy<Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class ProbabilisticDamagePolicy<Dim<3> >;
#endif
}