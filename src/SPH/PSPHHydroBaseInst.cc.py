text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "SPH/PSPHHydroBase.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class PSPHHydroBase< Dim< %(ndim)s > >;
}
"""
