text = """
//------------------------------------------------------------------------------
// Instantiations.
//------------------------------------------------------------------------------
#include "Damage/JohnsonCookDamage.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class PhysicsSpace::JohnsonCookDamage<Dim< %(ndim)s > >;
}
"""
