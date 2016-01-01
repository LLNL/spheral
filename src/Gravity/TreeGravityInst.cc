//------------------------------------------------------------------------------
// Explict instantiations.
//------------------------------------------------------------------------------
#include "TreeGravity.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace GravitySpace {

#ifdef SPHERAL2D
    template class TreeGravity<Dim<2> >;
#endif
#ifdef SPHERAL3D
    template class TreeGravity<Dim<3> >;
#endif

  }
}
