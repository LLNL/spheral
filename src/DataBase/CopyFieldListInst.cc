//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "CopyFieldList.cc"

namespace Spheral {
  template class CopyFieldList<Dim<1>, Dim<1>::Scalar>;
  template class CopyFieldList<Dim<1>, Dim<1>::Vector>;
  template class CopyFieldList<Dim<1>, Dim<1>::Vector3d>;
  template class CopyFieldList<Dim<1>, Dim<1>::Tensor>;
  template class CopyFieldList<Dim<1>, Dim<1>::SymTensor>;
                 
  template class CopyFieldList<Dim<2>, Dim<2>::Scalar>;
  template class CopyFieldList<Dim<2>, Dim<2>::Vector>;
  template class CopyFieldList<Dim<2>, Dim<2>::Vector3d>;
  template class CopyFieldList<Dim<2>, Dim<2>::Tensor>;
  template class CopyFieldList<Dim<2>, Dim<2>::SymTensor>;
                 
  template class CopyFieldList<Dim<3>, Dim<3>::Scalar>;
  template class CopyFieldList<Dim<3>, Dim<3>::Vector>;
  template class CopyFieldList<Dim<3>, Dim<3>::Vector3d>;
  template class CopyFieldList<Dim<3>, Dim<3>::Tensor>;
  template class CopyFieldList<Dim<3>, Dim<3>::SymTensor>;
}
