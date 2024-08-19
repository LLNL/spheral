//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Material/HelmholtzEquationOfState.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class HelmholtzEquationOfState< Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class HelmholtzEquationOfState< Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class HelmholtzEquationOfState< Dim<3> >;
#endif
}