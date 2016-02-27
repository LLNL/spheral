text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "RigidBoundary.cc"

namespace Spheral {
  namespace BoundarySpace {
    template class RigidBoundary< Dim< %(ndim)s > >;
  }
}
"""
