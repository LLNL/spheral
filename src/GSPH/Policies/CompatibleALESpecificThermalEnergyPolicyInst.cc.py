text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "Hydro/CompatibleMFVSpecificThermalEnergyPolicy.cc"

namespace Spheral {
  template class CompatibleMFVSpecificThermalEnergyPolicy<Dim< %(ndim)s > >;
}
"""
