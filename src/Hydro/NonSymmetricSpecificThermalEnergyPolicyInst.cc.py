text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"

#ifdef _OPENMP
#include "Hydro/NonSymmetricSpecificThermalEnergyPolicy_OpenMP.cc"
#else
#include "Hydro/NonSymmetricSpecificThermalEnergyPolicy.cc"
#endif

namespace Spheral {
  template class NonSymmetricSpecificThermalEnergyPolicy<Dim< %(ndim)s > >;
}
"""
