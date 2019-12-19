text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "NodeList/FluidNodeList.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class FluidNodeList< Dim< %(ndim)s > >;
}
"""
