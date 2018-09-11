text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Boundary/RigidBoundary.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace BoundarySpace {
    template class RigidBoundary< Dim< %(ndim)s > >;
  }
}
"""
