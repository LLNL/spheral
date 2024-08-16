//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Boundary/ConstantYVelocityBoundary.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_2D)
  template class ConstantYVelocityBoundary< Dim< 2 > >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class ConstantYVelocityBoundary< Dim< 3 > >;
#endif
}
