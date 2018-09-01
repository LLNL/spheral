text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "TreeDistributedBoundary.cc"

namespace Spheral {
  template class TreeDistributedBoundary< Dim< %(ndim)s > >;
}
"""
