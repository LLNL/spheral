//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "RK/RKCorrections.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
template class RKCorrections<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class RKCorrections<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class RKCorrections<Dim<3>>;
#endif
}