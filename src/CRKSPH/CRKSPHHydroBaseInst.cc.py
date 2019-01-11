text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "CRKSPHHydroBase.cc"
#include "Geometry/Dimension.hh"

#ifdef _OPENMP
#include "CRKSPHEvaluateDerivatives_OpenMP.cc"
#else
#include "CRKSPHEvaluateDerivatives.cc"
#endif

namespace Spheral {
template class CRKSPHHydroBase< Dim< %(ndim)s > >;
}
"""
