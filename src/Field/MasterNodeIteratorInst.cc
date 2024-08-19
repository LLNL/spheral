//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Field/MasterNodeIterator.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class MasterNodeIterator< Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class MasterNodeIterator< Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class MasterNodeIterator< Dim<3> >;
#endif
}