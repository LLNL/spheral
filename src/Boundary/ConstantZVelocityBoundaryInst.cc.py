text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "ConstantZVelocityBoundary.cc"

namespace Spheral {
  namespace BoundarySpace {
    template class ConstantZVelocityBoundary< Dim< %(ndim)s > >;
  }
}
"""
