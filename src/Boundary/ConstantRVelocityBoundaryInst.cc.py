text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "ConstantRVelocityBoundary.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class ConstantRVelocityBoundary< Dim< %(ndim)s > >;
}
"""
