text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "CRKSPH/CRKSPHHydroBase.cc"
#include "CRKSPH/CRKSPHEvaluateDerivatives.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
template class CRKSPHHydroBase< Dim< %(ndim)s > >;
}
"""
