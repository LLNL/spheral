text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "CRKSPHHydroBase.cc"

namespace Spheral {
  namespace CRKSPHSpace {
    template class CRKSPHHydroBase< Dim< %(ndim)s > >;
  }
}
"""
