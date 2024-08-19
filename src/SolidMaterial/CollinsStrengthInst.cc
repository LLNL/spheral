//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "SolidMaterial/CollinsStrength.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class CollinsStrength<Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class CollinsStrength<Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class CollinsStrength<Dim<3> >;
#endif
}