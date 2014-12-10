//------------------------------------------------------------------------------
// Explict instantiations.
//------------------------------------------------------------------------------
#include "TreeGravity.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace GravitySpace {

    template class TreeGravity<Dim<2> >;
    template class TreeGravity<Dim<3> >;

  }
}
