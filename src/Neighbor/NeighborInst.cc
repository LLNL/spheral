//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Neighbor/Neighbor.cc"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class Neighbor< Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class Neighbor< Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class Neighbor< Dim<3> >;
#endif
}