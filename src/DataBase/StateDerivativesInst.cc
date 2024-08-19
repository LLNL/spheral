//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "DataBase/StateDerivatives.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class StateDerivatives<Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class StateDerivatives<Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class StateDerivatives<Dim<3> >;
#endif
}