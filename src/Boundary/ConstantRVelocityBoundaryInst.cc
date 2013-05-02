//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "ConstantRVelocityBoundary.cc"

namespace Spheral {
  namespace BoundarySpace {
    template class ConstantRVelocityBoundary< Dim<1> >;
    template class ConstantRVelocityBoundary< Dim<2> >;
    template class ConstantRVelocityBoundary< Dim<3> >;
  }
}
