text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "Hydro/CompatibleDifferenceSpecificThermalEnergyPolicy.cc"

namespace Spheral {
  template class CompatibleDifferenceSpecificThermalEnergyPolicy<Dim< %(ndim)s > >;
}
"""
