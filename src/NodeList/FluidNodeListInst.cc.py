text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "FluidNodeList.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace NodeSpace {
    template class FluidNodeList< Dim< %(ndim)s > >;
  }
}
"""
