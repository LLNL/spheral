text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "SPH/SolidSPHHydroBase.cc"
#include "SPH/SolidSPHEvaluateDerivatives.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class SolidSPHHydroBase< Dim< %(ndim)s > >;
}
"""
