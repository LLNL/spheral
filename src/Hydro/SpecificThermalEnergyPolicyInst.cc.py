text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "SpecificThermalEnergyPolicy.cc"

namespace Spheral {
  template class SpecificThermalEnergyPolicy<Dim< %(ndim)s > >;
}
"""
