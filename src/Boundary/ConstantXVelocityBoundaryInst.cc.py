text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "ConstantXVelocityBoundary.cc"

namespace Spheral {
  namespace BoundarySpace {
    template class ConstantXVelocityBoundary< Dim< %(ndim)s > >;
  }
}
"""
