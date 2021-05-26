text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "RSPH/FSISpecificThermalEnergyPolicy.cc"

namespace Spheral {
  template class FSISpecificThermalEnergyPolicy<Dim< %(ndim)s > >;
}
"""
