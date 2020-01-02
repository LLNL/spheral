text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Hydro/NonSymmetricSpecificThermalEnergyPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class NonSymmetricSpecificThermalEnergyPolicy<Dim< %(ndim)s > >;
}
"""
