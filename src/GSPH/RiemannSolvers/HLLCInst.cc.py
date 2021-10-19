text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "GSPH/RiemannSolvers/HLLC.cc"

namespace Spheral {
  template class HLLC<Dim< %(ndim)s > >;
}
"""
