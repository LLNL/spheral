//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Kernel/PiGaussianKernel.hh"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class PiGaussianKernel< Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class PiGaussianKernel< Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class PiGaussianKernel< Dim<3> >;
#endif
}
