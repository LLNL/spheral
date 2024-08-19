//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "GSPH/MFVHydroBase.cc"
#include "GSPH/MFVEvaluateDerivatives.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class MFVHydroBase< Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class MFVHydroBase< Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class MFVHydroBase< Dim<3> >;
#endif
}