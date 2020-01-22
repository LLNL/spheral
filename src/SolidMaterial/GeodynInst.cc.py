text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SolidMaterial/Geodyn.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class Geodyn<Dim< %(ndim)s > >;
}
"""
