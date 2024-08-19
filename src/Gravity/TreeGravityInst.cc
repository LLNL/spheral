//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Gravity/TreeGravity.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_2D)
  template class TreeGravity<Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class TreeGravity<Dim<3> >;
#endif
}
