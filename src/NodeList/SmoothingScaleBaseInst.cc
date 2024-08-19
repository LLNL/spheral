//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "NodeList/SmoothingScaleBase.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class SmoothingScaleBase< Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class SmoothingScaleBase< Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class SmoothingScaleBase< Dim<3> >;
#endif
}