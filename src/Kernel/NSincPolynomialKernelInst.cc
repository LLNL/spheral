//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Kernel/NSincPolynomialKernel.cc"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class NSincPolynomialKernel< Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class NSincPolynomialKernel< Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class NSincPolynomialKernel< Dim<3> >;
#endif
}