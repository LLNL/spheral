text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "ConstantXVelocityBoundary.cc"

namespace Spheral {
  template class ConstantXVelocityBoundary< Dim< %(ndim)s > >;
}
"""
