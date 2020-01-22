text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "SPH/SPHHydroBase.cc"
#include "SPH/SPHEvaluateDerivatives.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class SPHHydroBase< Dim< %(ndim)s > >;
}
"""
