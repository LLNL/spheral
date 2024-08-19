//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Kernel/WendlandC6Kernel.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class WendlandC6Kernel< Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class WendlandC6Kernel< Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class WendlandC6Kernel< Dim<3> >;
#endif
}
