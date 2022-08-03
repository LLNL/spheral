text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "SPH/SPHHydroBase.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class SPHHydroBase<%(Dim)s>;
}
"""
