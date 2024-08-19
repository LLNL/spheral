//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "FlatConnectivity.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class FlatConnectivity<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class FlatConnectivity<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class FlatConnectivity<Dim<3>>;
#endif
}