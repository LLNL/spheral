text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "CRKSPHHydroBase.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace CRKSPHSpace {
    template class CRKSPHHydroBase< Dim< %(ndim)s > >;
  }
}
"""
