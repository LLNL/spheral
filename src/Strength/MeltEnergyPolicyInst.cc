//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Strength/MeltEnergyPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class MeltEnergyPolicy<Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class MeltEnergyPolicy<Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class MeltEnergyPolicy<Dim<3> >;
#endif
}