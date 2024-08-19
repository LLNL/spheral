//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "SPH/SPHHydroBase.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class SPHHydroBase<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class SPHHydroBase<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class SPHHydroBase<Dim<3>>;
#endif
}