//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "DEM/SolidBoundary/DEMBoundaryPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class DEMBoundaryPolicy<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class DEMBoundaryPolicy<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class DEMBoundaryPolicy<Dim<3>>;
#endif
}