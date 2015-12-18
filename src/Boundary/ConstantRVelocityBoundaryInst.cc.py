text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "ConstantRVelocityBoundary.cc"

namespace Spheral {
  namespace BoundarySpace {
    template class ConstantRVelocityBoundary< Dim< %(ndim)s > >;
  }
}
"""
