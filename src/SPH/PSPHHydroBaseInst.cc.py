text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "PSPHHydroBase.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace SPHSpace {
    template class PSPHHydroBase< Dim< %(ndim)s > >;
  }
}
"""
