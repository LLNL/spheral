//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Porosity/PorosityModel.cc"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class PorosityModel<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class PorosityModel<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class PorosityModel<Dim<3>>;
#endif
}