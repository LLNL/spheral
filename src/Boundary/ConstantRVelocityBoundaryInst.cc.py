text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Boundary/ConstantRVelocityBoundary.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class ConstantRVelocityBoundary< Dim< %(ndim)s > >;
}
"""
