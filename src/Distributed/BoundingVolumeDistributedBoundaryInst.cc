//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "BoundingVolumeDistributedBoundary.cc"

namespace Spheral {
  namespace BoundarySpace {
    template class BoundingVolumeDistributedBoundary< Dim<1> >;
    template class BoundingVolumeDistributedBoundary< Dim<2> >;
    template class BoundingVolumeDistributedBoundary< Dim<3> >;
  }
}
