text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Gravity/CompatibleGravitationalVelocityPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class CompatibleGravitationalVelocityPolicy<Dim< %(ndim)s > >;
}
"""
