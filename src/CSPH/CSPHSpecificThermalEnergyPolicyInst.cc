//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "CSPHSpecificThermalEnergyPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class CSPHSpecificThermalEnergyPolicy<Dim<1> >;
  template class CSPHSpecificThermalEnergyPolicy<Dim<2> >;
  template class CSPHSpecificThermalEnergyPolicy<Dim<3> >;
}
