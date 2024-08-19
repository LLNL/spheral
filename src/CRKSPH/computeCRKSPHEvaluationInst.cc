//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "CRKSPH/computeCRKSPHEvaluation.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
template void computeCRKSPHEvaluation(const ConnectivityMap<Dim<1> >&,
                                     const TableKernel<Dim<1> >&,
                                     const FieldList<Dim<1>, Dim<1>::Scalar>&,
                                     const FieldList<Dim<1>, Dim<1>::Vector>&,
                                     const FieldList<Dim<1>, Dim<1>::SymTensor>&,
                                     size_t, const int, Dim<1>::Vector,
                                     const bool, Dim<1>::Scalar&, Dim<1>::Vector&);
#endif

#if defined(SPHERAL_ENABLE_2D)
template void computeCRKSPHEvaluation(const ConnectivityMap<Dim<2> >&,
                                     const TableKernel<Dim<2> >&,
                                     const FieldList<Dim<2>, Dim<2>::Scalar>&,
                                     const FieldList<Dim<2>, Dim<2>::Vector>&,
                                     const FieldList<Dim<2>, Dim<2>::SymTensor>&,
                                     size_t, const int, Dim<2>::Vector,
                                     const bool, Dim<2>::Scalar&, Dim<2>::Vector&);
#endif

#if defined(SPHERAL_ENABLE_3D)
template void computeCRKSPHEvaluation(const ConnectivityMap<Dim<3> >&,
                                     const TableKernel<Dim<3> >&,
                                     const FieldList<Dim<3>, Dim<3>::Scalar>&,
                                     const FieldList<Dim<3>, Dim<3>::Vector>&,
                                     const FieldList<Dim<3>, Dim<3>::SymTensor>&,
                                     size_t, const int, Dim<3>::Vector,
                                     const bool, Dim<3>::Scalar&, Dim<3>::Vector&);
#endif
}