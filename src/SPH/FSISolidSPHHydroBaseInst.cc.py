text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "SPH/FSISolidSPHHydroBase.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class FSISolidSPHHydroBase< Dim< %(ndim)s > >;
}
"""
