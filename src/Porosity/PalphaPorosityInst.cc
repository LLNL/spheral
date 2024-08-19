//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Porosity/PalphaPorosity.cc"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class PalphaPorosity<Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class PalphaPorosity<Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class PalphaPorosity<Dim<3> >;
#endif
}