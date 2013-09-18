//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "DistributedBoundary.cc"

namespace Spheral {
  namespace BoundarySpace {
    template class DistributedBoundary< Dim<1> >;
    template class DistributedBoundary< Dim<2> >;
    template class DistributedBoundary< Dim<3> >;
  }
}
