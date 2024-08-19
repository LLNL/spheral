//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "NodeList/NodeListRegistrar.cc"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class NodeListRegistrar<Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class NodeListRegistrar<Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class NodeListRegistrar<Dim<3> >;
#endif
}
