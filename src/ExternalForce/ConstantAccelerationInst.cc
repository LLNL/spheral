//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "ConstantAcceleration.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace PhysicsSpace {
    template class ConstantAcceleration< Dim<1> >;
    template class ConstantAcceleration< Dim<2> >;
    template class ConstantAcceleration< Dim<3> >;
  }
}
