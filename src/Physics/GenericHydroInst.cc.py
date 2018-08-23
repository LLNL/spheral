text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "Physics/GenericHydro.cc"

namespace Spheral {
  namespace PhysicsSpace {
    template class GenericHydro< Dim< %(ndim)s > >;
  }
}
"""
