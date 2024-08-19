//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Distributed/DistributedBoundary.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class DistributedBoundary< Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class DistributedBoundary< Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class DistributedBoundary< Dim<3> >;
#endif
}