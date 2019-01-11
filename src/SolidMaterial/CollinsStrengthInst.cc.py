text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "CollinsStrength.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class CollinsStrength<Dim< %(ndim)s > >;
}
"""
