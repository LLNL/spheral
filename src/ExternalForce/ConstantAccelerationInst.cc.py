text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "ExternalForce/ConstantAcceleration.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace PhysicsSpace {
    template class ConstantAcceleration< Dim< %(ndim)s > >;
  }
}
"""
