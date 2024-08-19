//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Material/PolytropicEquationOfState.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class PolytropicEquationOfState< Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class PolytropicEquationOfState< Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class PolytropicEquationOfState< Dim<3> >;
#endif
}