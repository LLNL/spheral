//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Material/IsothermalEquationOfState.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class IsothermalEquationOfState< Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class IsothermalEquationOfState< Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class IsothermalEquationOfState< Dim<3> >;
#endif
}