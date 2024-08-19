//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Kernel/WendlandC4Kernel.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class WendlandC4Kernel< Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class WendlandC4Kernel< Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class WendlandC4Kernel< Dim<3> >;
#endif
}
