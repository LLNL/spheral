//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "KernelIntegrator.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class KernelIntegrator<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class KernelIntegrator<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class KernelIntegrator<Dim<3>>;
#endif
}