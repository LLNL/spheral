text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "FluidNodeList.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class FluidNodeList< Dim< %(ndim)s > >;
}
"""
