//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "SVPH/computeSVPHCorrections.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)

template
void
computeSVPHCorrections(const ConnectivityMap<Dim<1> >& connectivityMap,
                       const TableKernel<Dim<1> >& W,
                       const FieldList<Dim<1>, Dim<1>::Scalar>& volume,
                       const FieldList<Dim<1>, Dim<1>::Vector>& position,
                       const FieldList<Dim<1>, Dim<1>::SymTensor>& H,
                       Field<Dim<1>, Dim<1>::Scalar>& A,
                       Field<Dim<1>, Dim<1>::Vector>& B,
                       Field<Dim<1>, Dim<1>::Tensor>& gradB);

template
void
computeSVPHCorrections(const ConnectivityMap<Dim<1> >& connectivityMap,
                       const TableKernel<Dim<1> >& W,
                       const FieldList<Dim<1>, Dim<1>::Scalar>& volume,
                       const FieldList<Dim<1>, Dim<1>::Vector>& position,
                       const FieldList<Dim<1>, Dim<1>::SymTensor>& H,
                       FieldList<Dim<1>, Dim<1>::Scalar>& A,
                       FieldList<Dim<1>, Dim<1>::Vector>& B,
                       FieldList<Dim<1>, Dim<1>::Tensor>& gradB);

#endif

#if defined(SPHERAL_ENABLE_2D)

template
void
computeSVPHCorrections(const ConnectivityMap<Dim<2> >& connectivityMap,
                       const TableKernel<Dim<2> >& W,
                       const FieldList<Dim<2>, Dim<2>::Scalar>& volume,
                       const FieldList<Dim<2>, Dim<2>::Vector>& position,
                       const FieldList<Dim<2>, Dim<2>::SymTensor>& H,
                       Field<Dim<2>, Dim<2>::Scalar>& A,
                       Field<Dim<2>, Dim<2>::Vector>& B,
                       Field<Dim<2>, Dim<2>::Tensor>& gradB);

template
void
computeSVPHCorrections(const ConnectivityMap<Dim<2> >& connectivityMap,
                       const TableKernel<Dim<2> >& W,
                       const FieldList<Dim<2>, Dim<2>::Scalar>& volume,
                       const FieldList<Dim<2>, Dim<2>::Vector>& position,
                       const FieldList<Dim<2>, Dim<2>::SymTensor>& H,
                       FieldList<Dim<2>, Dim<2>::Scalar>& A,
                       FieldList<Dim<2>, Dim<2>::Vector>& B,
                       FieldList<Dim<2>, Dim<2>::Tensor>& gradB);

#endif

#if defined(SPHERAL_ENABLE_3D)

template
void
computeSVPHCorrections(const ConnectivityMap<Dim<3> >& connectivityMap,
                       const TableKernel<Dim<3> >& W,
                       const FieldList<Dim<3>, Dim<3>::Scalar>& volume,
                       const FieldList<Dim<3>, Dim<3>::Vector>& position,
                       const FieldList<Dim<3>, Dim<3>::SymTensor>& H,
                       Field<Dim<3>, Dim<3>::Scalar>& A,
                       Field<Dim<3>, Dim<3>::Vector>& B,
                       Field<Dim<3>, Dim<3>::Tensor>& gradB);

template
void
computeSVPHCorrections(const ConnectivityMap<Dim<3> >& connectivityMap,
                       const TableKernel<Dim<3> >& W,
                       const FieldList<Dim<3>, Dim<3>::Scalar>& volume,
                       const FieldList<Dim<3>, Dim<3>::Vector>& position,
                       const FieldList<Dim<3>, Dim<3>::SymTensor>& H,
                       FieldList<Dim<3>, Dim<3>::Scalar>& A,
                       FieldList<Dim<3>, Dim<3>::Vector>& B,
                       FieldList<Dim<3>, Dim<3>::Tensor>& gradB);

#endif
}