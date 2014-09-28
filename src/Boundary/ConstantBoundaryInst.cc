//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "ConstantBoundary.cc"

namespace Spheral {
  namespace BoundarySpace {
    template class ConstantBoundary< Dim<1> >;
    template class ConstantBoundary< Dim<2> >;
    template class ConstantBoundary< Dim<3> >;
  }
}
