//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Gravity/Tree.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_2D)
  template class Tree<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class Tree<Dim<3>>;
#endif
}
