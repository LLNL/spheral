//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "NodeList/NodeList.cc"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class NodeList< Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class NodeList< Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class NodeList< Dim<3> >;
#endif
}