//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "NestedGridDistributedBoundary.cc"

namespace Spheral {
  namespace BoundarySpace {
    template class NestedGridDistributedBoundary< Dim<1> >;
    template class NestedGridDistributedBoundary< Dim<2> >;
    template class NestedGridDistributedBoundary< Dim<3> >;
  }
}
