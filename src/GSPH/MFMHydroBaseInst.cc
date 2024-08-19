//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "GSPH/MFMHydroBase.cc"
#include "GSPH/MFMEvaluateDerivatives.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class MFMHydroBase< Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class MFMHydroBase< Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class MFMHydroBase< Dim<3> >;
#endif
}