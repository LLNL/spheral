text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "SPHHydroBase.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace SPHSpace {
    template class SPHHydroBase< Dim< %(ndim)s > >;
  }
}
"""
