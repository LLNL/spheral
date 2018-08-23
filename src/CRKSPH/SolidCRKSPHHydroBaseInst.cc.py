text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "CRKSPH/SolidCRKSPHHydroBase.cc"

#ifdef _OPENMP
#include "CRKSPH/SolidCRKSPHEvaluateDerivatives_OpenMP.cc"
#else
#include "CRKSPH/SolidCRKSPHEvaluateDerivatives.cc"
#endif

namespace Spheral {
  namespace CRKSPHSpace {
    template class SolidCRKSPHHydroBase< Dim< %(ndim)s > >;
  }
}
"""
