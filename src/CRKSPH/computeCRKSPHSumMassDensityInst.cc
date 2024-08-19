//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Boundary/Boundary.hh"
#include "CRKSPH/computeCRKSPHSumMassDensity.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
template void computeCRKSPHSumMassDensity(const ConnectivityMap<Dim<1> >&,
                                          const TableKernel<Dim<1> >&,
                                          const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                          const FieldList<Dim<1>, Dim<1>::Scalar>& mass,
                                          const FieldList<Dim<1>, Dim<1>::Scalar>& vol,
                                          const FieldList<Dim<1>, Dim<1>::SymTensor>& H,
                                          FieldList<Dim<1>, Dim<1>::Scalar>& massDensity);
#endif

#if defined(SPHERAL_ENABLE_2D)
template void computeCRKSPHSumMassDensity(const ConnectivityMap<Dim<2> >&,
                                          const TableKernel<Dim<2> >&,
                                          const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                          const FieldList<Dim<2>, Dim<2>::Scalar>& mass,
                                          const FieldList<Dim<2>, Dim<2>::Scalar>& vol,
                                          const FieldList<Dim<2>, Dim<2>::SymTensor>& H,
                                          FieldList<Dim<2>, Dim<2>::Scalar>& massDensity);
#endif

#if defined(SPHERAL_ENABLE_3D)
template void computeCRKSPHSumMassDensity(const ConnectivityMap<Dim<3> >&,
                                          const TableKernel<Dim<3> >&,
                                          const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                          const FieldList<Dim<3>, Dim<3>::Scalar>& mass,
                                          const FieldList<Dim<3>, Dim<3>::Scalar>& vol,
                                          const FieldList<Dim<3>, Dim<3>::SymTensor>& H,
                                          FieldList<Dim<3>, Dim<3>::Scalar>& massDensity);
#endif
}