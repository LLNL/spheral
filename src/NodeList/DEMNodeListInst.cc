//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "NodeList/DEMNodeList.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class DEMNodeList< Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class DEMNodeList< Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class DEMNodeList< Dim<3> >;
#endif
}