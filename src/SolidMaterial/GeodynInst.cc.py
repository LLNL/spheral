text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SolidMaterial/Geodyn.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace SolidMaterial {
    template class Geodyn<Dim< %(ndim)s > >;
  }
}
"""
