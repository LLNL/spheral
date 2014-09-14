//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "GenericHydro.cc"

namespace Spheral {
  namespace PhysicsSpace {
    template class GenericHydro< Dim<1> >;
    template class GenericHydro< Dim<2> >;
    template class GenericHydro< Dim<3> >;
  }
}
