//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "RK/ReproducingKernelMethods.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
template class ReproducingKernelMethods<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class ReproducingKernelMethods<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class ReproducingKernelMethods<Dim<3>>;
#endif
}