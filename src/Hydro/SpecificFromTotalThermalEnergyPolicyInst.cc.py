text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Hydro/SpecificFromTotalThermalEnergyPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class SpecificFromTotalThermalEnergyPolicy<Dim< %(ndim)s > >;
}
"""
