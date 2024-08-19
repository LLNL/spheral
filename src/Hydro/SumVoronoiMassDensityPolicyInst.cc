//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Hydro/SumVoronoiMassDensityPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class SumVoronoiMassDensityPolicy<Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class SumVoronoiMassDensityPolicy<Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class SumVoronoiMassDensityPolicy<Dim<3> >;
#endif
}