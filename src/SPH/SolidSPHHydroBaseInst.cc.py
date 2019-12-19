text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "SolidSPHHydroBase.cc"
#include "SolidSPHEvaluateDerivatives.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class SolidSPHHydroBase< Dim< %(ndim)s > >;
}
"""
