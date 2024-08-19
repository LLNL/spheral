//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Kernel/Kernel.hh"
#include "Kernel/WendlandC2Kernel.hh"
#include <math.h>

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class WendlandC2Kernel< Dim<1>  >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class WendlandC2Kernel< Dim<2>  >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class WendlandC2Kernel< Dim<3>  >;
#endif
}