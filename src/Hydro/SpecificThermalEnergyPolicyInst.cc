//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SpecificThermalEnergyPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class SpecificThermalEnergyPolicy<Dim<1> >;
  template class SpecificThermalEnergyPolicy<Dim<2> >;
  template class SpecificThermalEnergyPolicy<Dim<3> >;
}
