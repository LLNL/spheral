//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Damage/TensorDamageModel.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class TensorDamageModel<Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class TensorDamageModel<Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class TensorDamageModel<Dim<3> >;
#endif
}