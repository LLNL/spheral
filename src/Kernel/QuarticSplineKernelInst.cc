//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Kernel/QuarticSplineKernel.hh"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class QuarticSplineKernel< Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class QuarticSplineKernel< Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class QuarticSplineKernel< Dim<3> >;
#endif
}
