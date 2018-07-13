text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "MeltEnergyPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class MeltEnergyPolicy<Dim< %(ndim)s > >;
}
"""
