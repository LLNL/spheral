//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Kernel/Kernel.hh"
#include "Kernel/HatKernel.hh"
#include <math.h>

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class HatKernel< Dim<1>  >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class HatKernel< Dim<2>  >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class HatKernel< Dim<3>  >;
#endif
}