text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "Hydro/SpecificThermalEnergyPolicy.cc"

namespace Spheral {
  template class SpecificThermalEnergyPolicy<Dim< %(ndim)s > >;
}
"""
