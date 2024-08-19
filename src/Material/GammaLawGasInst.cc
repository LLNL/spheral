//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Material/GammaLawGas.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class GammaLawGas< Dim<1>  >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class GammaLawGas< Dim<2>  >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class GammaLawGas< Dim<3>  >;
#endif
}