text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "NodeList/DEMNodeList.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class DEMNodeList< Dim< %(ndim)s > >;
}
"""
