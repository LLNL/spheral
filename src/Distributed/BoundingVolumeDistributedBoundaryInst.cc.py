text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "BoundingVolumeDistributedBoundary.cc"

namespace Spheral {
  namespace BoundarySpace {
    template class BoundingVolumeDistributedBoundary< Dim< %(ndim)s > >;
  }
}
"""
