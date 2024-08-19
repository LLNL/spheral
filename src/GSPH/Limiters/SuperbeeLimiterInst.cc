//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "GSPH/Limiters/SuperbeeLimiter.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class SuperbeeLimiter<Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class SuperbeeLimiter<Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class SuperbeeLimiter<Dim<3> >;
#endif
}