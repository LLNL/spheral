//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "GSPH/GSPHHydroBase.cc"
#include "GSPH/GSPHEvaluateDerivatives.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class GSPHHydroBase< Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class GSPHHydroBase< Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class GSPHHydroBase< Dim<3> >;
#endif
}