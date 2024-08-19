//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Hydro/SecondMomentHourglassControl.cc"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class SecondMomentHourglassControl< Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class SecondMomentHourglassControl< Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class SecondMomentHourglassControl< Dim<3> >;
#endif
}