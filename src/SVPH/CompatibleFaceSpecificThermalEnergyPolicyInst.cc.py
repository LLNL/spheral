text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "CompatibleFaceSpecificThermalEnergyPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class CompatibleFaceSpecificThermalEnergyPolicy<Dim< %(ndim)s > >;
}
"""
