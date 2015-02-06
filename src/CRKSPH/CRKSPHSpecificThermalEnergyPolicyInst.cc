//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "CRKSPHSpecificThermalEnergyPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class CRKSPHSpecificThermalEnergyPolicy<Dim<1> >;
  template class CRKSPHSpecificThermalEnergyPolicy<Dim<2> >;
  template class CRKSPHSpecificThermalEnergyPolicy<Dim<3> >;
}
