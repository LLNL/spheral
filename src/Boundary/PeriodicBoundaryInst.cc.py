text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "Boundary/PeriodicBoundary.cc"

namespace Spheral {
  template class PeriodicBoundary< Dim< %(ndim)s > >;
}
"""
