//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Hydro/SpecificFromTotalThermalEnergyPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class SpecificFromTotalThermalEnergyPolicy<Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class SpecificFromTotalThermalEnergyPolicy<Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class SpecificFromTotalThermalEnergyPolicy<Dim<3> >;
#endif
}