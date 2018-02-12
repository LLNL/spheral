text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "CompatibleGravitationalVelocityPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class CompatibleGravitationalVelocityPolicy<Dim< %(ndim)s > >;
}
"""
