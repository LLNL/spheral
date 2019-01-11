text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "CRKSPH/CRKSPHHydroBase.cc"
#include "Geometry/Dimension.hh"

#ifdef _OPENMP
#include "CRKSPH/CRKSPHEvaluateDerivatives_OpenMP.cc"
#else
#include "CRKSPH/CRKSPHEvaluateDerivatives.cc"
#endif

namespace Spheral {
template class CRKSPHHydroBase< Dim< %(ndim)s > >;
}
"""
