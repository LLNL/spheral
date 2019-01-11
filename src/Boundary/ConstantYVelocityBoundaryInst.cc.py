text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "ConstantYVelocityBoundary.cc"

namespace Spheral {
  template class ConstantYVelocityBoundary< Dim< %(ndim)s > >;
}
"""
