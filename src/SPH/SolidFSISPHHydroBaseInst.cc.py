text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "SPH/SolidFSISPHHydroBase.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class SolidFSISPHHydroBase< Dim< %(ndim)s > >;
}
"""
