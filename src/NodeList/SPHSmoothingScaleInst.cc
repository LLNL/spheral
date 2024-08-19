//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "NodeList/SPHSmoothingScale.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class SPHSmoothingScale<Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class SPHSmoothingScale<Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class SPHSmoothingScale<Dim<3> >;
#endif
}
