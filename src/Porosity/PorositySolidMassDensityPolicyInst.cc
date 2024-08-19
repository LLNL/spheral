//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Porosity/PorositySolidMassDensityPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class PorositySolidMassDensityPolicy<Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class PorositySolidMassDensityPolicy<Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class PorositySolidMassDensityPolicy<Dim<3> >;
#endif
}