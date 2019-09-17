text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "SPHHydroBase.cc"
#include "SPHEvaluateDerivatives.cc"
#include "Geometry/Dimension.hh"

// #ifdef _OPENMP
// #include "SPHEvaluateDerivatives_OpenMP.cc"
// #else
// #include "SPHEvaluateDerivatives.cc"
// #endif

namespace Spheral {
  template class SPHHydroBase< Dim< %(ndim)s > >;
}
"""
