text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SolidMaterial/JohnsonCookStrength.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class JohnsonCookStrength<Dim< %(ndim)s > >;
}
"""
