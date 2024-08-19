//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "RK/computeHVolumes.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template void computeHVolumes(const Dim<1>::Scalar nPerh,
                                const FieldList<Dim<1>, Dim<1>::SymTensor>&,
                                FieldList<Dim<1>, Dim<1>::Scalar>&);
#endif

#if defined(SPHERAL_ENABLE_2D)
  template void computeHVolumes(const Dim<2>::Scalar nPerh,
                                const FieldList<Dim<2>, Dim<2>::SymTensor>&,
                                FieldList<Dim<2>, Dim<2>::Scalar>&);
#endif

#if defined(SPHERAL_ENABLE_3D)
  template void computeHVolumes(const Dim<3>::Scalar nPerh,
                                const FieldList<Dim<3>, Dim<3>::SymTensor>&,
                                FieldList<Dim<3>, Dim<3>::Scalar>&);
#endif
}