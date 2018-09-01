text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "SolidCRKSPHHydroBase.cc"

#ifdef _OPENMP
#include "SolidCRKSPHEvaluateDerivatives_OpenMP.cc"
#else
#include "SolidCRKSPHEvaluateDerivatives.cc"
#endif

namespace Spheral {
template class SolidCRKSPHHydroBase< Dim< %(ndim)s > >;
}
"""
