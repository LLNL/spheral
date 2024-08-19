//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Strength/ShearModulusPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class ShearModulusPolicy<Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class ShearModulusPolicy<Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class ShearModulusPolicy<Dim<3> >;
#endif
}