text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "Boundary/ConstantYVelocityBoundary.cc"

namespace Spheral {
  template class ConstantYVelocityBoundary< Dim< %(ndim)s > >;
}
"""
