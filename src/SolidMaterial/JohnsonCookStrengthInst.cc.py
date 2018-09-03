text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "JohnsonCookStrength.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class JohnsonCookStrength<Dim< %(ndim)s > >;
}
"""
