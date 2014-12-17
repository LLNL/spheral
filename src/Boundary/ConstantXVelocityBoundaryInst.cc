//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "ConstantXVelocityBoundary.cc"

namespace Spheral {
  namespace BoundarySpace {
    template class ConstantXVelocityBoundary< Dim<1> >;
    template class ConstantXVelocityBoundary< Dim<2> >;
    template class ConstantXVelocityBoundary< Dim<3> >;
  }
}
