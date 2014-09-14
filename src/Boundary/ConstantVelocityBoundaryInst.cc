//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "ConstantVelocityBoundary.cc"

namespace Spheral {
  namespace BoundarySpace {
    template class ConstantVelocityBoundary< Dim<1> >;
    template class ConstantVelocityBoundary< Dim<2> >;
    template class ConstantVelocityBoundary< Dim<3> >;
  }
}
