text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "NestedGridDistributedBoundary.cc"

namespace Spheral {
  namespace BoundarySpace {
    template class NestedGridDistributedBoundary< Dim< %(ndim)s > >;
  }
}
"""
