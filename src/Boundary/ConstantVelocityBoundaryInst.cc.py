text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "ConstantVelocityBoundary.cc"

namespace Spheral {
  namespace BoundarySpace {
    template class ConstantVelocityBoundary< Dim< %(ndim)s > >;
  }
}
"""
