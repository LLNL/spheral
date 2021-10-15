text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "GSPH/GSPHSpecificThermalEnergyPolicy.cc"

namespace Spheral {
  template class GSPHSpecificThermalEnergyPolicy<Dim< %(ndim)s > >;
}
"""
