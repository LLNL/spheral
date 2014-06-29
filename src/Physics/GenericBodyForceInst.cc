//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "GenericBodyForce.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace PhysicsSpace {
    template class GenericBodyForce< Dim<1> >;
    template class GenericBodyForce< Dim<2> >;
    template class GenericBodyForce< Dim<3> >;
  }
}
