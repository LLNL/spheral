text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SVPH/CompatibleFaceSpecificThermalEnergyPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class CompatibleFaceSpecificThermalEnergyPolicy<Dim< %(ndim)s > >;
}
"""
