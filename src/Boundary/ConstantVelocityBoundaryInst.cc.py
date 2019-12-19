text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Boundary/ConstantVelocityBoundary.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class ConstantVelocityBoundary< Dim< %(ndim)s > >;
}
"""
