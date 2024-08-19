//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Damage/TensorStrainPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class Spheral::TensorStrainPolicy<Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class Spheral::TensorStrainPolicy<Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class Spheral::TensorStrainPolicy<Dim<3> >;
#endif
}