text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "GSPH/RiemannSolvers/GHLLC.cc"

namespace Spheral {
  template class GHLLC<Dim< %(ndim)s > >;
}
"""
