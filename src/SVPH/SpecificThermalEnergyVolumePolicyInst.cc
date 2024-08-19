//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "SVPH/SpecificThermalEnergyVolumePolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class SpecificThermalEnergyVolumePolicy<Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class SpecificThermalEnergyVolumePolicy<Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class SpecificThermalEnergyVolumePolicy<Dim<3> >;
#endif
}