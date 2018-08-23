text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SVPH/SpecificThermalEnergyVolumePolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class SpecificThermalEnergyVolumePolicy<Dim< %(ndim)s > >;
}
"""
