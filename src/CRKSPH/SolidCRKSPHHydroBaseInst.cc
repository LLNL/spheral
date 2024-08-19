//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "CRKSPH/SolidCRKSPHHydroBase.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
template class SolidCRKSPHHydroBase< Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class SolidCRKSPHHydroBase< Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class SolidCRKSPHHydroBase< Dim<3> >;
#endif
}