//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "SPHGravity/SPHGravity.cc"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class SPHGravity<Dim<1> >;
} // end namespace Spheral

"""
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class SPHGravity<Dim<2> >;
} // end namespace Spheral

"""
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class SPHGravity<Dim<3> >;
} // end namespace Spheral

"""
#endif
}