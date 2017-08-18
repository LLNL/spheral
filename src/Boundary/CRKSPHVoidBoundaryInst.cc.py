text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "CRKSPHVoidBoundary.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace BoundarySpace {
    template class CRKSPHVoidBoundary< Dim< %(ndim)s > >;
  }
}
"""
