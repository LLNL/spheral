text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "SPH/SPHHydroBase.cc"
#include "Geometry/Dimension.hh"

#ifdef _OPENMP
#include "SPH/SPHEvaluateDerivatives_OpenMP.cc"
#else
#include "SPH/SPHEvaluateDerivatives.cc"
#endif

namespace Spheral {
  template class SPHHydroBase< Dim< %(ndim)s > >;
}
"""
