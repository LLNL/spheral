text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"

#ifdef _OPENMP
#include "NonSymmetricSpecificThermalEnergyPolicy_OpenMP.cc"
#else
#include "NonSymmetricSpecificThermalEnergyPolicy.cc"
#endif

namespace Spheral {
  template class NonSymmetricSpecificThermalEnergyPolicy<Dim< %(ndim)s > >;
}
"""
