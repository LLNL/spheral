//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "BilinearIndex.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class BilinearIndex<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class BilinearIndex<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class BilinearIndex<Dim<3>>;
#endif
}