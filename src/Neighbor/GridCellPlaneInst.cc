//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Neighbor/GridCellPlane.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class GridCellPlane< Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class GridCellPlane< Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class GridCellPlane< Dim<3> >;
#endif
}