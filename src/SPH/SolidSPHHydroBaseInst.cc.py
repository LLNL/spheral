text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "SolidSPHHydroBase.cc"
#include "Geometry/Dimension.hh"

#ifdef _OPENMP
#include "SolidSPHEvaluateDerivatives_OpenMP.cc"
#else
#include "SolidSPHEvaluateDerivatives.cc"
#endif

namespace Spheral {
  namespace SPHSpace {
    template class SolidSPHHydroBase< Dim< %(ndim)s > >;
  }
}
"""
