text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "SPH/SPHHydroBase.cc"
#include "Geometry/Dimension.hh"

#ifdef _OPENMP
#include "SPHEvaluateDerivatives_OpenMP.cc"
#else
#include "SPHEvaluateDerivatives.cc"
#endif

namespace Spheral {
  template class SPHHydroBase< Dim< %(ndim)s > >;
}
"""
