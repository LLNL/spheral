//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "SolidMaterial/StrengthModel.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class StrengthModel<Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class StrengthModel<Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class StrengthModel<Dim<3> >;
#endif
}