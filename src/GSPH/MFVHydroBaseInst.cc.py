text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "GSPH/MFVHydroBase.cc"
#include "GSPH/MFVEvaluateDerivatives.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class MFVHydroBase< Dim< %(ndim)s > >;
}
"""
