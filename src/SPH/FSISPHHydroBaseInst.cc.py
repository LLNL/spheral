text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "SPH/FSISPHHydroBase.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class FSISPHHydroBase< Dim< %(ndim)s > >;
}
"""
