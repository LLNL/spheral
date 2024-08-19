//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Gravity/NBodyGravity.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)

  template class NBodyGravity<Dim<1> >;

} // end namespace Spheral

"""
#endif

#if defined(SPHERAL_ENABLE_2D)

  template class NBodyGravity<Dim<2> >;

} // end namespace Spheral

"""
#endif

#if defined(SPHERAL_ENABLE_3D)

  template class NBodyGravity<Dim<3> >;

} // end namespace Spheral

"""
#endif
}