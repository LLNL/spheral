//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "NodeList/ASPHSmoothingScale.cc"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class ASPHSmoothingScale<Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class ASPHSmoothingScale<Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class ASPHSmoothingScale<Dim<3> >;
#endif
}