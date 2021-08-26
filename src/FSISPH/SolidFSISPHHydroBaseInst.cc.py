text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "FSISPH/SolidFSISPHHydroBase.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class SolidFSISPHHydroBase< Dim< %(ndim)s > >;
}
"""
