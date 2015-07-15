//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "ReplaceBoundedFieldList.cc"

namespace Spheral {
  template class ReplaceBoundedFieldList<Dim<1>, Dim<1>::Scalar>;
  template class ReplaceBoundedFieldList<Dim<1>, Dim<1>::Vector>;
  template class ReplaceBoundedFieldList<Dim<1>, Dim<1>::Vector3d>;
  template class ReplaceBoundedFieldList<Dim<1>, Dim<1>::Tensor>;
  template class ReplaceBoundedFieldList<Dim<1>, Dim<1>::SymTensor>;
  template class ReplaceBoundedFieldList<Dim<1>, Dim<1>::SymTensor, Dim<1>::Scalar>;
                 
  template class ReplaceBoundedFieldList<Dim<2>, Dim<2>::Scalar>;
  template class ReplaceBoundedFieldList<Dim<2>, Dim<2>::Vector>;
  template class ReplaceBoundedFieldList<Dim<2>, Dim<2>::Vector3d>;
  template class ReplaceBoundedFieldList<Dim<2>, Dim<2>::Tensor>;
  template class ReplaceBoundedFieldList<Dim<2>, Dim<2>::SymTensor>;
  template class ReplaceBoundedFieldList<Dim<2>, Dim<2>::SymTensor, Dim<2>::Scalar>;
                 
  template class ReplaceBoundedFieldList<Dim<3>, Dim<3>::Scalar>;
  template class ReplaceBoundedFieldList<Dim<3>, Dim<3>::Vector>;
  template class ReplaceBoundedFieldList<Dim<3>, Dim<3>::Vector3d>;
  template class ReplaceBoundedFieldList<Dim<3>, Dim<3>::Tensor>;
  template class ReplaceBoundedFieldList<Dim<3>, Dim<3>::SymTensor>;
  template class ReplaceBoundedFieldList<Dim<3>, Dim<3>::SymTensor, Dim<3>::Scalar>;
}
