text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "CRKSPHHydroBase.cc"
#include "CRKSPHEvaluateDerivatives.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
template class CRKSPHHydroBase< Dim< %(ndim)s > >;
}
"""
