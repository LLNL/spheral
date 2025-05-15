text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "OverlapNodeMap.cc"

namespace Spheral {
  template class OverlapNodeMap<Dim<%(ndim)s>>;
}
"""
