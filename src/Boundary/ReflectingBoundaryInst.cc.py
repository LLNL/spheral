text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "ReflectingBoundary.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class ReflectingBoundary< Dim< %(ndim)s > >;
}
"""
