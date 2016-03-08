text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "ConstantBoundary.cc"

namespace Spheral {
  namespace BoundarySpace {
    template class ConstantBoundary< Dim< %(ndim)s > >;
  }
}
"""
