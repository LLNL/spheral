//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Physics/GenericHydro.cc"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class GenericHydro< Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class GenericHydro< Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class GenericHydro< Dim<3> >;
#endif
}