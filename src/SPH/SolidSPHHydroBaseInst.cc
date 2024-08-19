//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "SPH/SolidSPHHydroBase.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class SolidSPHHydroBase< Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class SolidSPHHydroBase< Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class SolidSPHHydroBase< Dim<3> >;
#endif
}