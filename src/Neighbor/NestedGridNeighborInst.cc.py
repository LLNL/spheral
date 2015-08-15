text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "NestedGridNeighbor.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace NeighborSpace {
    template class NestedGridNeighbor< Dim< %(ndim)s > >;
  }
}
"""
