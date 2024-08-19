//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Distributed/PeanoHilbertOrderRedistributeNodes.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class PeanoHilbertOrderRedistributeNodes< Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class PeanoHilbertOrderRedistributeNodes< Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class PeanoHilbertOrderRedistributeNodes< Dim<3> >;
#endif
}