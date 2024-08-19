//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "NodeList/FixedSmoothingScale.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class FixedSmoothingScale<Dim< 1 > >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class FixedSmoothingScale<Dim< 2 > >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class FixedSmoothingScale<Dim< 3 > >;
#endif
}
