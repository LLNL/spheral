text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "ReflectingBoundary.cc"

namespace Spheral {
  namespace BoundarySpace {
    template class ReflectingBoundary< Dim< %(ndim)s > >;
  }
}
"""
