text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "NodeList/NeighborNodeList.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class NeighborNodeList< Dim< %(ndim)s > >;
}
"""
