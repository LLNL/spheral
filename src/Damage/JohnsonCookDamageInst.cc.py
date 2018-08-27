text = """
//------------------------------------------------------------------------------
// Instantiations.
//------------------------------------------------------------------------------
#include "JohnsonCookDamage.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class PhysicsSpace::JohnsonCookDamage<Dim< %(ndim)s > >;
}
"""
