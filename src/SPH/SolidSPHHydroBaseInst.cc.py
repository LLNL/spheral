text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "SPH/SolidSPHHydroBase.cc"
#include "Geometry/Dimension.hh"

#ifdef _OPENMP
#include "SPH/SolidSPHEvaluateDerivatives_OpenMP.cc"
#else
#include "SPH/SolidSPHEvaluateDerivatives.cc"
#endif

namespace Spheral {
  namespace SPHSpace {
    template class SolidSPHHydroBase< Dim< %(ndim)s > >;
  }
}
"""
