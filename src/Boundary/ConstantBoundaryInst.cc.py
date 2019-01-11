text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Boundary/ConstantBoundary.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class ConstantBoundary< Dim< %(ndim)s > >;
}
"""
