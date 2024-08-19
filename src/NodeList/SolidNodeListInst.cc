//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "NodeList/SolidNodeList.cc"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class SolidNodeList< Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class SolidNodeList< Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class SolidNodeList< Dim<3> >;
#endif
}