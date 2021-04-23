text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "FSISPH/FSISPHHydroBase.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class FSISPHHydroBase< Dim< %(ndim)s > >;
}
"""
