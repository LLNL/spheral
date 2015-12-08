text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "GenericBodyForce.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace PhysicsSpace {
    template class GenericBodyForce< Dim< %(ndim)s > >;
  }
}
"""
