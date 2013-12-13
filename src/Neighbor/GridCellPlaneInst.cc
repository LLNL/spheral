//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "GridCellPlane.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace NeighborSpace {
    template class GridCellPlane< Dim<1> >;
    template class GridCellPlane< Dim<2> >;
    template class GridCellPlane< Dim<3> >;
  }
}
