//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "CRKSPH/CRKSPHHydroBase.cc"
#include "CRKSPH/CRKSPHEvaluateDerivatives.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
template class CRKSPHHydroBase< Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class CRKSPHHydroBase< Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class CRKSPHHydroBase< Dim<3> >;
#endif
}