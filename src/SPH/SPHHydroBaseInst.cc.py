text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "SPHHydroBase.cc"
#include "SPHEvaluateDerivatives.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class SPHHydroBase< Dim< %(ndim)s > >;
}
"""
