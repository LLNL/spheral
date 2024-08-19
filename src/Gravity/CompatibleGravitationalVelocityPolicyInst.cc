//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Gravity/CompatibleGravitationalVelocityPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class CompatibleGravitationalVelocityPolicy<Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class CompatibleGravitationalVelocityPolicy<Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class CompatibleGravitationalVelocityPolicy<Dim<3> >;
#endif
}