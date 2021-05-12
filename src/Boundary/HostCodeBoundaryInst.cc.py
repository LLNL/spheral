text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Boundary/HostCodeBoundary.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class HostCodeBoundary< Dim< %(ndim)s > >;
}
"""
