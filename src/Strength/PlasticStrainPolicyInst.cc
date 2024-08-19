//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Strength/PlasticStrainPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class PlasticStrainPolicy<Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class PlasticStrainPolicy<Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class PlasticStrainPolicy<Dim<3> >;
#endif
}