text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "SolidCRKSPHHydroBase.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace CRKSPHSpace {
    template class SolidCRKSPHHydroBase< Dim< %(ndim)s > >;
  }
}
"""
