text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "NodeMap.cc"

namespace Spheral {
  template class NodeMap<Dim<%(ndim)s>>;
}
"""
