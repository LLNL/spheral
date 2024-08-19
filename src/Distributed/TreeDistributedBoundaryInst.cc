//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Distributed/TreeDistributedBoundary.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class TreeDistributedBoundary< Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class TreeDistributedBoundary< Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class TreeDistributedBoundary< Dim<3> >;
#endif
}