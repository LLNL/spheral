text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SolidMaterial/CollinsStrength.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class CollinsStrength<Dim< %(ndim)s > >;
}
"""
