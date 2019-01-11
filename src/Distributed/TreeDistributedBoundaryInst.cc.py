text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "Distributed/TreeDistributedBoundary.cc"

namespace Spheral {
  template class TreeDistributedBoundary< Dim< %(ndim)s > >;
}
"""
