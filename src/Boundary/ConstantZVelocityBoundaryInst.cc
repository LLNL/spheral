//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "ConstantZVelocityBoundary.cc"

namespace Spheral {
  template class ConstantZVelocityBoundary< Dim<3> >;
}
