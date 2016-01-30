text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "PeriodicBoundary.cc"

namespace Spheral {
  namespace BoundarySpace {
    template class PeriodicBoundary< Dim< %(ndim)s > >;
  }
}
"""
