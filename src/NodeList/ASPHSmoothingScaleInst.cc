//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"

// Define CPP macros for specializations in the .cc file.
#if defined(SPHERAL_ENABLE_1D)
#define SPHERAL1DINSTANTIATION
#endif

#if defined(SPHERAL_ENABLE_2D)
#define SPHERAL2DINSTANTIATION
#endif

#if defined(SPHERAL_ENABLE_3D)
#define SPHERAL3DINSTANTIATION
#endif

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
