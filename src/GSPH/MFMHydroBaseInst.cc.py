text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "GSPH/MFMHydroBase.cc"
#include "GSPH/MFMEvaluateDerivatives.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class MFMHydroBase< Dim< %(ndim)s > >;
}
"""
