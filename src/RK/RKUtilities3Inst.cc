//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "RK/RKUtilities.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
template class RKUtilities<Dim<1>, RKOrder::CubicOrder>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class RKUtilities<Dim<2>, RKOrder::CubicOrder>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class RKUtilities<Dim<3>, RKOrder::CubicOrder>;
#endif
}
