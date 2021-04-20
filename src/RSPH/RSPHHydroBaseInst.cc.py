text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "RSPH/RSPHHydroBase.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class RSPHHydroBase< Dim< %(ndim)s > >;
}
"""
