//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "GSPH/RiemannSolvers/RiemannSolverBase.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class RiemannSolverBase<Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class RiemannSolverBase<Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class RiemannSolverBase<Dim<3> >;
#endif
}