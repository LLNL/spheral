//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "ReflectingBoundary.cc"

namespace Spheral {
  namespace BoundarySpace {
    template class ReflectingBoundary< Dim<1> >;
    template class ReflectingBoundary< Dim<2> >;
    template class ReflectingBoundary< Dim<3> >;
  }
}
