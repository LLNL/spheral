text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "Distributed/TreeDistributedBoundary.cc"

namespace Spheral {
  namespace BoundarySpace {
    template class TreeDistributedBoundary< Dim< %(ndim)s > >;
  }
}
"""
