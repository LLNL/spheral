//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "RK/computeHullVolumes.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template void computeHullVolumes(const ConnectivityMap<Dim<1> >&,
                                   const Dim<1>::Scalar kernelExtent,
                                   const FieldList<Dim<1>, Dim<1>::Vector>&,
                                   const FieldList<Dim<1>, Dim<1>::SymTensor>&,
                                   FieldList<Dim<1>, Dim<1>::Scalar>&);
#endif

#if defined(SPHERAL_ENABLE_2D)
  template void computeHullVolumes(const ConnectivityMap<Dim<2> >&,
                                   const Dim<2>::Scalar kernelExtent,
                                   const FieldList<Dim<2>, Dim<2>::Vector>&,
                                   const FieldList<Dim<2>, Dim<2>::SymTensor>&,
                                   FieldList<Dim<2>, Dim<2>::Scalar>&);
#endif

#if defined(SPHERAL_ENABLE_3D)
  template void computeHullVolumes(const ConnectivityMap<Dim<3> >&,
                                   const Dim<3>::Scalar kernelExtent,
                                   const FieldList<Dim<3>, Dim<3>::Vector>&,
                                   const FieldList<Dim<3>, Dim<3>::SymTensor>&,
                                   FieldList<Dim<3>, Dim<3>::Scalar>&);
#endif
}