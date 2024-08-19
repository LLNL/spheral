//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "SolidMaterial/iSALEROCKStrength.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class iSALEROCKStrength<Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class iSALEROCKStrength<Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class iSALEROCKStrength<Dim<3> >;
#endif
}