text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SolidMaterial/ANEOS.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace SolidMaterial {
    template class ANEOS<Dim< %(ndim)s > >;
  }
}
"""
