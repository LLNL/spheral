//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Kernel/SincKernel.hh"
#include <math.h>

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class SincKernel< Dim<1>  >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class SincKernel< Dim<2>  >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class SincKernel< Dim<3>  >;
#endif
}