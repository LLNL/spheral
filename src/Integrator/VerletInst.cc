//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Integrator/Verlet.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class Verlet< Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class Verlet< Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class Verlet< Dim<3> >;
#endif
}