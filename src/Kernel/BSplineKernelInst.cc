//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Kernel/Kernel.hh"
#include "Kernel/BSplineKernel.hh"
#include <math.h>

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class BSplineKernel<Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class BSplineKernel<Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class BSplineKernel<Dim<3> >;
#endif
}