//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Neighbor/TreeNeighbor.cc"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class TreeNeighbor< Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class TreeNeighbor< Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class TreeNeighbor< Dim<3> >;
#endif
}