//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Utilities/computeShepardsInterpolation.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
template
FieldList<Dim<1>, Dim<1>::Scalar>
computeShepardsInterpolation(const FieldList<Dim<1>, Dim<1>::Scalar>& fieldList,
                             const ConnectivityMap<Dim<1> >&,
                             const TableKernel<Dim<1> >&,
                             const FieldList<Dim<1>, Dim<1>::Vector>& position,
                             const FieldList<Dim<1>, Dim<1>::SymTensor>& H,
                             const FieldList<Dim<1>, Dim<1>::Scalar>& weight);

template
FieldList<Dim<1>, Dim<1>::Vector>
computeShepardsInterpolation(const FieldList<Dim<1>, Dim<1>::Vector>& fieldList,
                             const ConnectivityMap<Dim<1> >&,
                             const TableKernel<Dim<1> >&,
                             const FieldList<Dim<1>, Dim<1>::Vector>& position,
                             const FieldList<Dim<1>, Dim<1>::SymTensor>& H,
                             const FieldList<Dim<1>, Dim<1>::Scalar>& weight);

template
FieldList<Dim<1>, Dim<1>::Tensor>
computeShepardsInterpolation(const FieldList<Dim<1>, Dim<1>::Tensor>& fieldList,
                             const ConnectivityMap<Dim<1> >&,
                             const TableKernel<Dim<1> >&,
                             const FieldList<Dim<1>, Dim<1>::Vector>& position,
                             const FieldList<Dim<1>, Dim<1>::SymTensor>& H,
                             const FieldList<Dim<1>, Dim<1>::Scalar>& weight);

template
FieldList<Dim<1>, Dim<1>::SymTensor>
computeShepardsInterpolation(const FieldList<Dim<1>, Dim<1>::SymTensor>& fieldList,
                             const ConnectivityMap<Dim<1> >&,
                             const TableKernel<Dim<1> >&,
                             const FieldList<Dim<1>, Dim<1>::Vector>& position,
                             const FieldList<Dim<1>, Dim<1>::SymTensor>& H,
                             const FieldList<Dim<1>, Dim<1>::Scalar>& weight);
#endif

#if defined(SPHERAL_ENABLE_2D)
template
FieldList<Dim<2>, Dim<2>::Scalar>
computeShepardsInterpolation(const FieldList<Dim<2>, Dim<2>::Scalar>& fieldList,
                             const ConnectivityMap<Dim<2> >&,
                             const TableKernel<Dim<2> >&,
                             const FieldList<Dim<2>, Dim<2>::Vector>& position,
                             const FieldList<Dim<2>, Dim<2>::SymTensor>& H,
                             const FieldList<Dim<2>, Dim<2>::Scalar>& weight);

template
FieldList<Dim<2>, Dim<2>::Vector>
computeShepardsInterpolation(const FieldList<Dim<2>, Dim<2>::Vector>& fieldList,
                             const ConnectivityMap<Dim<2> >&,
                             const TableKernel<Dim<2> >&,
                             const FieldList<Dim<2>, Dim<2>::Vector>& position,
                             const FieldList<Dim<2>, Dim<2>::SymTensor>& H,
                             const FieldList<Dim<2>, Dim<2>::Scalar>& weight);

template
FieldList<Dim<2>, Dim<2>::Tensor>
computeShepardsInterpolation(const FieldList<Dim<2>, Dim<2>::Tensor>& fieldList,
                             const ConnectivityMap<Dim<2> >&,
                             const TableKernel<Dim<2> >&,
                             const FieldList<Dim<2>, Dim<2>::Vector>& position,
                             const FieldList<Dim<2>, Dim<2>::SymTensor>& H,
                             const FieldList<Dim<2>, Dim<2>::Scalar>& weight);

template
FieldList<Dim<2>, Dim<2>::SymTensor>
computeShepardsInterpolation(const FieldList<Dim<2>, Dim<2>::SymTensor>& fieldList,
                             const ConnectivityMap<Dim<2> >&,
                             const TableKernel<Dim<2> >&,
                             const FieldList<Dim<2>, Dim<2>::Vector>& position,
                             const FieldList<Dim<2>, Dim<2>::SymTensor>& H,
                             const FieldList<Dim<2>, Dim<2>::Scalar>& weight);
#endif

#if defined(SPHERAL_ENABLE_3D)
template
FieldList<Dim<3>, Dim<3>::Scalar>
computeShepardsInterpolation(const FieldList<Dim<3>, Dim<3>::Scalar>& fieldList,
                             const ConnectivityMap<Dim<3> >&,
                             const TableKernel<Dim<3> >&,
                             const FieldList<Dim<3>, Dim<3>::Vector>& position,
                             const FieldList<Dim<3>, Dim<3>::SymTensor>& H,
                             const FieldList<Dim<3>, Dim<3>::Scalar>& weight);

template
FieldList<Dim<3>, Dim<3>::Vector>
computeShepardsInterpolation(const FieldList<Dim<3>, Dim<3>::Vector>& fieldList,
                             const ConnectivityMap<Dim<3> >&,
                             const TableKernel<Dim<3> >&,
                             const FieldList<Dim<3>, Dim<3>::Vector>& position,
                             const FieldList<Dim<3>, Dim<3>::SymTensor>& H,
                             const FieldList<Dim<3>, Dim<3>::Scalar>& weight);

template
FieldList<Dim<3>, Dim<3>::Tensor>
computeShepardsInterpolation(const FieldList<Dim<3>, Dim<3>::Tensor>& fieldList,
                             const ConnectivityMap<Dim<3> >&,
                             const TableKernel<Dim<3> >&,
                             const FieldList<Dim<3>, Dim<3>::Vector>& position,
                             const FieldList<Dim<3>, Dim<3>::SymTensor>& H,
                             const FieldList<Dim<3>, Dim<3>::Scalar>& weight);

template
FieldList<Dim<3>, Dim<3>::SymTensor>
computeShepardsInterpolation(const FieldList<Dim<3>, Dim<3>::SymTensor>& fieldList,
                             const ConnectivityMap<Dim<3> >&,
                             const TableKernel<Dim<3> >&,
                             const FieldList<Dim<3>, Dim<3>::Vector>& position,
                             const FieldList<Dim<3>, Dim<3>::SymTensor>& H,
                             const FieldList<Dim<3>, Dim<3>::Scalar>& weight);
#endif
}
