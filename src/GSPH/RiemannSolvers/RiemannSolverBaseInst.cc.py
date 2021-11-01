text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "GSPH/RiemannSolvers/RiemannSolverBase.cc"

namespace Spheral {
  template class RiemannSolverBase<Dim< %(ndim)s > >;
}
"""
