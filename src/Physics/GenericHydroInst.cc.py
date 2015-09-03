text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "GenericHydro.cc"

namespace Spheral {
  namespace PhysicsSpace {
    template class GenericHydro< Dim< %(ndim)s > >;
  }
}
"""
