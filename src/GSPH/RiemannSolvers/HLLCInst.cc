//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "GSPH/RiemannSolvers/HLLC.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class HLLC<Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class HLLC<Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class HLLC<Dim<3> >;
#endif
}