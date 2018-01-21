text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "LinearAcceleration.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace PhysicsSpace {
    template class LinearAcceleration< Dim< %(ndim)s > >;
  }
}
"""
