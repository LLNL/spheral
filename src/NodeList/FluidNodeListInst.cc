//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "NodeList/FluidNodeList.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class FluidNodeList< Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class FluidNodeList< Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class FluidNodeList< Dim<3> >;
#endif
}