//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "ExternalForce/PointPotential.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class PointPotential< Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class PointPotential< Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class PointPotential< Dim<3> >;
#endif
}