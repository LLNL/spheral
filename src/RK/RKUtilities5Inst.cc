//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "RK/RKUtilities.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
template class RKUtilities<Dim<1>, """
    text += """RKOrder::%(order)s>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class RKUtilities<Dim<2>, """
    text += """RKOrder::%(order)s>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class RKUtilities<Dim<3>, """
    text += """RKOrder::%(order)s>;
#endif
}