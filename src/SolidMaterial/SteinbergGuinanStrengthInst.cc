//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "SolidMaterial/SteinbergGuinanStrength.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class SteinbergGuinanStrength<Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class SteinbergGuinanStrength<Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class SteinbergGuinanStrength<Dim<3> >;
#endif
}