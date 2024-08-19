//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Utilities/overlayRemapFields.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
template void
overlayRemapFields(const std::vector<Boundary<Dim<1>>*>& boundaries,
                   const std::vector<Field<Dim<1>, Dim<1>::Scalar>*>& scalarDonorFields,
                   const std::vector<Field<Dim<1>, Dim<1>::Vector>*>& vectorDonorFields,
                   const std::vector<Field<Dim<1>, Dim<1>::Tensor>*>& tensorDonorFields,
                   const std::vector<Field<Dim<1>, Dim<1>::SymTensor>*>& symTensorDonorFields,
                   std::vector<Field<Dim<1>, Dim<1>::Scalar>*>& scalarAcceptorFields,
                   std::vector<Field<Dim<1>, Dim<1>::Vector>*>& vectorAcceptorFields,
                   std::vector<Field<Dim<1>, Dim<1>::Tensor>*>& tensorAcceptorFields,
                   std::vector<Field<Dim<1>, Dim<1>::SymTensor>*>& symTensorAcceptorFields);
#endif

#if defined(SPHERAL_ENABLE_2D)
template void
overlayRemapFields(const std::vector<Boundary<Dim<2>>*>& boundaries,
                   const std::vector<Field<Dim<2>, Dim<2>::Scalar>*>& scalarDonorFields,
                   const std::vector<Field<Dim<2>, Dim<2>::Vector>*>& vectorDonorFields,
                   const std::vector<Field<Dim<2>, Dim<2>::Tensor>*>& tensorDonorFields,
                   const std::vector<Field<Dim<2>, Dim<2>::SymTensor>*>& symTensorDonorFields,
                   std::vector<Field<Dim<2>, Dim<2>::Scalar>*>& scalarAcceptorFields,
                   std::vector<Field<Dim<2>, Dim<2>::Vector>*>& vectorAcceptorFields,
                   std::vector<Field<Dim<2>, Dim<2>::Tensor>*>& tensorAcceptorFields,
                   std::vector<Field<Dim<2>, Dim<2>::SymTensor>*>& symTensorAcceptorFields);
#endif

#if defined(SPHERAL_ENABLE_3D)
template void
overlayRemapFields(const std::vector<Boundary<Dim<3>>*>& boundaries,
                   const std::vector<Field<Dim<3>, Dim<3>::Scalar>*>& scalarDonorFields,
                   const std::vector<Field<Dim<3>, Dim<3>::Vector>*>& vectorDonorFields,
                   const std::vector<Field<Dim<3>, Dim<3>::Tensor>*>& tensorDonorFields,
                   const std::vector<Field<Dim<3>, Dim<3>::SymTensor>*>& symTensorDonorFields,
                   std::vector<Field<Dim<3>, Dim<3>::Scalar>*>& scalarAcceptorFields,
                   std::vector<Field<Dim<3>, Dim<3>::Vector>*>& vectorAcceptorFields,
                   std::vector<Field<Dim<3>, Dim<3>::Tensor>*>& tensorAcceptorFields,
                   std::vector<Field<Dim<3>, Dim<3>::SymTensor>*>& symTensorAcceptorFields);
#endif
}