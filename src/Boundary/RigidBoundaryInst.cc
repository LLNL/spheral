//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "RigidBoundary.cc"

namespace Spheral {
  namespace BoundarySpace {
    template class RigidBoundary< Dim<1> >;
    template class RigidBoundary< Dim<2> >;
    template class RigidBoundary< Dim<3> >;
  }
}
