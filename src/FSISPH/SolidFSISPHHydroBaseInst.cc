//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "FSISPH/SolidFSISPHHydroBase.cc"
#include "FSISPH/SolidFSISPHEvaluateDerivatives.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class SolidFSISPHHydroBase< Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class SolidFSISPHHydroBase< Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class SolidFSISPHHydroBase< Dim<3> >;
#endif
}