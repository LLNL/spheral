text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "TaylorSPH/TaylorSPHHydroBase.cc"

namespace Spheral {
  template class TaylorSPHHydroBase< Dim< %(ndim)s > >;
}
"""
