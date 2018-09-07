text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "CRKSPHVoidBoundary.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class CRKSPHVoidBoundary< Dim< %(ndim)s > >;
}
"""
