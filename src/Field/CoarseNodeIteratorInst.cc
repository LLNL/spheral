//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Field/CoarseNodeIterator.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class CoarseNodeIterator< Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class CoarseNodeIterator< Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class CoarseNodeIterator< Dim<3> >;
#endif
}