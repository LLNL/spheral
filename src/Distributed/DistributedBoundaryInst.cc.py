text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "Distributed/DistributedBoundary.cc"

namespace Spheral {
  template class DistributedBoundary< Dim< %(ndim)s > >;
}
"""
