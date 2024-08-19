//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Hydro/ThirdMomentHourglassControl.cc"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class ThirdMomentHourglassControl< Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class ThirdMomentHourglassControl< Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class ThirdMomentHourglassControl< Dim<3> >;
#endif
}