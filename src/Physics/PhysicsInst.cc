//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "Physics.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace PhysicsSpace {
    template class Physics<Dim<1> >;
    template class Physics<Dim<2> >;
    template class Physics<Dim<3> >;
  }
}
