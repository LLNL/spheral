//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "GSPH/Limiters/LimiterBase.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class LimiterBase<Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class LimiterBase<Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class LimiterBase<Dim<3> >;
#endif
}