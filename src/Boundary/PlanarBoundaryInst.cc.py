text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "PlanarBoundary.cc"

namespace Spheral {
  namespace BoundarySpace {
    template class PlanarBoundary< Dim< %(ndim)s > >;
  }
}
"""
