//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "CRKSPH/CRKSPHVariant.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
template class CRKSPHVariant< Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class CRKSPHVariant< Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class CRKSPHVariant< Dim<3> >;
#endif
}