//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Utilities/medianPosition.cc"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template Dim<1>::Vector medianPosition(vector<Dim<1>::Vector>& positions);
#endif

#if defined(SPHERAL_ENABLE_2D)
  template Dim<2>::Vector medianPosition(vector<Dim<2>::Vector>& positions);
#endif

#if defined(SPHERAL_ENABLE_3D)
  template Dim<3>::Vector medianPosition(vector<Dim<3>::Vector>& positions);
#endif
}