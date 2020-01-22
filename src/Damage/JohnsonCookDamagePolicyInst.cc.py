text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "Damage/JohnsonCookDamagePolicy.cc"

namespace Spheral {
  template class JohnsonCookDamagePolicy<Dim< %(ndim)s > >;
}
"""
