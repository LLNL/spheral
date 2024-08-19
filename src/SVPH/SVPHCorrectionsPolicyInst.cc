//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "SVPH/SVPHCorrectionsPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class SVPHCorrectionsPolicy<Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class SVPHCorrectionsPolicy<Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class SVPHCorrectionsPolicy<Dim<3> >;
#endif
}