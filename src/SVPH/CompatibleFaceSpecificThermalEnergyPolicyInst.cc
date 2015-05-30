//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "CompatibleFaceSpecificThermalEnergyPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class CompatibleFaceSpecificThermalEnergyPolicy<Dim<1> >;
  template class CompatibleFaceSpecificThermalEnergyPolicy<Dim<2> >;
  template class CompatibleFaceSpecificThermalEnergyPolicy<Dim<3> >;
}
