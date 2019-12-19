text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Boundary/RigidBoundary.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class RigidBoundary< Dim< %(ndim)s > >;
}
"""
