//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Kernel/Kernel.hh"
#include "Kernel/W4SplineKernel.hh"
#include <math.h>

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class W4SplineKernel< Dim<1>  >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class W4SplineKernel< Dim<2>  >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class W4SplineKernel< Dim<3>  >;
#endif
}