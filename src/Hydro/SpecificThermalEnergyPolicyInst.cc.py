text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"

#ifdef _OPENMP
#include "Hydro/SpecificThermalEnergyPolicy_OpenMP.cc"
#else
#include "Hydro/SpecificThermalEnergyPolicy.cc"
#endif

namespace Spheral {
  template class SpecificThermalEnergyPolicy<Dim< %(ndim)s > >;
}
"""
