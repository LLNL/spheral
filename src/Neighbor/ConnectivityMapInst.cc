//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Neighbor/ConnectivityMap.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class ConnectivityMap<Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class ConnectivityMap<Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class ConnectivityMap<Dim<3> >;
#endif
}
