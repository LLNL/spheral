text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "TaylorSPHHydroBase.cc"

namespace Spheral {
  namespace TaylorSPHSpace {
    template class TaylorSPHHydroBase< Dim< %(ndim)s > >;
  }
}
"""
