//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "SolidMaterial/OsborneEquationOfState.cc"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class OsborneEquationOfState<Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class OsborneEquationOfState<Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class OsborneEquationOfState<Dim<3> >;
#endif
}