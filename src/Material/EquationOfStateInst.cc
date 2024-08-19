//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Material/EquationOfState.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class EquationOfState< Dim<1>  >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class EquationOfState< Dim<2>  >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class EquationOfState< Dim<3>  >;
#endif
}