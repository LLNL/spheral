text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "NestedGridDistributedBoundary.cc"

namespace Spheral {
  template class NestedGridDistributedBoundary< Dim< %(ndim)s > >;
}
"""
