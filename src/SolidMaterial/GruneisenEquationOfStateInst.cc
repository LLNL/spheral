//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "SolidMaterial/GruneisenEquationOfState.cc"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class GruneisenEquationOfState<Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class GruneisenEquationOfState<Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class GruneisenEquationOfState<Dim<3> >;
#endif
}