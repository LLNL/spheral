//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "DEM/LinearSpringDEM.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class LinearSpringDEM< Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class LinearSpringDEM< Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class LinearSpringDEM< Dim<3> >;
#endif
}