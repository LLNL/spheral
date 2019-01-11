text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "TaylorSPHHydroBase.cc"

namespace Spheral {
  template class TaylorSPHHydroBase< Dim< %(ndim)s > >;
}
"""
