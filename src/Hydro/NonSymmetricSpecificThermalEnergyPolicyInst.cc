//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "NonSymmetricSpecificThermalEnergyPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class NonSymmetricSpecificThermalEnergyPolicy<Dim<1> >;
  template class NonSymmetricSpecificThermalEnergyPolicy<Dim<2> >;
  template class NonSymmetricSpecificThermalEnergyPolicy<Dim<3> >;
}
