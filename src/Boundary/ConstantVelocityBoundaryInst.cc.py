text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "ConstantVelocityBoundary.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class ConstantVelocityBoundary< Dim< %(ndim)s > >;
}
"""
