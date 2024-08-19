//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Integrator/CheapSynchronousRK2.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class CheapSynchronousRK2< Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class CheapSynchronousRK2< Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class CheapSynchronousRK2< Dim<3> >;
#endif
}