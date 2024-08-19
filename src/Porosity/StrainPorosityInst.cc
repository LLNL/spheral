//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Porosity/StrainPorosity.cc"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class StrainPorosity<Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class StrainPorosity<Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class StrainPorosity<Dim<3> >;
#endif
}