//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Physics/Physics.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class Physics<Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class Physics<Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class Physics<Dim<3> >;
#endif
}