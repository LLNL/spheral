//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "GSPH/Policies/MFVIncrementSpecificThermalEnergyPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class MFVIncrementSpecificThermalEnergyPolicy<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class MFVIncrementSpecificThermalEnergyPolicy<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class MFVIncrementSpecificThermalEnergyPolicy<Dim<3>>;
#endif
}