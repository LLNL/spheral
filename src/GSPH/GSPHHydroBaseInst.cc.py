text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "GSPH/GSPHHydroBase.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class GSPHHydroBase< Dim< %(ndim)s > >;
}
"""
