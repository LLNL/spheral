text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "SolidSPHHydroBase.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace SPHSpace {
    template class SolidSPHHydroBase< Dim< %(ndim)s > >;
  }
}
"""
