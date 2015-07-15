//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "PlanarBoundary.cc"

namespace Spheral {
  namespace BoundarySpace {
    template class PlanarBoundary< Dim<1> >;
    template class PlanarBoundary< Dim<2> >;
    template class PlanarBoundary< Dim<3> >;
  }
}
