text = """
//------------------------------------------------------------------------------
// Instantiations.
//------------------------------------------------------------------------------
#include "JohnsonCookDamageBase.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class PhysicsSpace::JohnsonCookDamageBase<Dim< %(ndim)s > >;
}
"""
