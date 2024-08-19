//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Boundary/ConstantZVelocityBoundary.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_3D)
  template class ConstantZVelocityBoundary< Dim<3> >;
#endif
}
