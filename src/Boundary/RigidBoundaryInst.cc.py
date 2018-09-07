text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "RigidBoundary.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class RigidBoundary< Dim< %(ndim)s > >;
}
"""
