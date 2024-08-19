//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "FVCRKH/FVCRKHydroBase.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class FVCRKHydroBase< Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class FVCRKHydroBase< Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class FVCRKHydroBase< Dim<3> >;
#endif
}
