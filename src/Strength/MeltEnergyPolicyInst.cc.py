text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Strength/MeltEnergyPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class MeltEnergyPolicy<Dim< %(ndim)s > >;
}
"""
