//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Hydro/PressurePolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class PressurePolicy<Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class PressurePolicy<Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class PressurePolicy<Dim<3> >;
#endif
}