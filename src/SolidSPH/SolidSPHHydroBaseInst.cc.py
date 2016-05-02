text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "SolidSPHHydroBase.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace SolidSPHSpace {
    template class SolidSPHHydroBase< Dim< %(ndim)s > >;
  }
}
"""
