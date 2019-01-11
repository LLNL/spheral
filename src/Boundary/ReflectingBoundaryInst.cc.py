text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Boundary/ReflectingBoundary.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class ReflectingBoundary< Dim< %(ndim)s > >;
}
"""
