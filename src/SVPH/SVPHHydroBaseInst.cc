//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "SVPH/SVPHHydroBase.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class SVPHHydroBase< Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class SVPHHydroBase< Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class SVPHHydroBase< Dim<3> >;
#endif
}
