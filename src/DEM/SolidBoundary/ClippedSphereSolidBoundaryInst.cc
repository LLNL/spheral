//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "DEM/SolidBoundary/ClippedSphereSolidBoundary.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class ClippedSphereSolidBoundary< Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class ClippedSphereSolidBoundary< Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class ClippedSphereSolidBoundary< Dim<3> >;
#endif
}