//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "LinearAcceleration.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace PhysicsSpace {
    template class LinearAcceleration< Dim<1> >;
    template class LinearAcceleration< Dim<2> >;
    template class LinearAcceleration< Dim<3> >;
  }
}
