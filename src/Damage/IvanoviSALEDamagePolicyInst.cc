//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Damage/IvanoviSALEDamagePolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class IvanoviSALEDamagePolicy<Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class IvanoviSALEDamagePolicy<Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class IvanoviSALEDamagePolicy<Dim<3> >;
#endif
}