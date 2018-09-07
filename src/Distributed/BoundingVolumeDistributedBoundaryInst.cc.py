text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "BoundingVolumeDistributedBoundary.cc"

namespace Spheral {
  template class BoundingVolumeDistributedBoundary< Dim< %(ndim)s > >;
}
"""
