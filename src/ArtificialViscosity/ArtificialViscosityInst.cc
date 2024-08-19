//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "ArtificialViscosity/ArtificialViscosity.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class ArtificialViscosity< Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class ArtificialViscosity< Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class ArtificialViscosity< Dim<3> >;
#endif
}