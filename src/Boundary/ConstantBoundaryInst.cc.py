text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "ConstantBoundary.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class ConstantBoundary< Dim< %(ndim)s > >;
}
"""
