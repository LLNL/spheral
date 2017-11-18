text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "TreeDistributedBoundary.cc"

namespace Spheral {
  namespace BoundarySpace {
    template class TreeDistributedBoundary< Dim< %(ndim)s > >;
  }
}
"""
