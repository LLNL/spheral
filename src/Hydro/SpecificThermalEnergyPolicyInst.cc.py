text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"

#ifdef _OPENMP
#include "SpecificThermalEnergyPolicy_OpenMP.cc"
#else
#include "SpecificThermalEnergyPolicy.cc"
#endif

namespace Spheral {
  template class SpecificThermalEnergyPolicy<Dim< %(ndim)s > >;
}
"""
