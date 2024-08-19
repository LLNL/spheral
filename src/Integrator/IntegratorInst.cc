//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Integrator/Integrator.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class Integrator< Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class Integrator< Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class Integrator< Dim<3> >;
#endif
}