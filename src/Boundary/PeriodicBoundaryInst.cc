//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "PeriodicBoundary.cc"

namespace Spheral {
  namespace BoundarySpace {
    template class PeriodicBoundary< Dim<1> >;
    template class PeriodicBoundary< Dim<2> >;
    template class PeriodicBoundary< Dim<3> >;
  }
}
