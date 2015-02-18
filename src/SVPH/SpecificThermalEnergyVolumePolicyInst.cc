//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SpecificThermalEnergyVolumePolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class SpecificThermalEnergyVolumePolicy<Dim<1> >;
  template class SpecificThermalEnergyVolumePolicy<Dim<2> >;
  template class SpecificThermalEnergyVolumePolicy<Dim<3> >;
}
