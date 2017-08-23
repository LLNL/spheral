text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "PlanarBoundary.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace BoundarySpace {
    template class PlanarBoundary< Dim< %(ndim)s > >;
  }
}
"""
