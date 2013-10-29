//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "ReplaceFieldList.cc"

namespace Spheral {
  template class ReplaceFieldList<Dim<1>, Dim<1>::Scalar>;
  template class ReplaceFieldList<Dim<1>, Dim<1>::Vector>;
  template class ReplaceFieldList<Dim<1>, Dim<1>::Vector3d>;
  template class ReplaceFieldList<Dim<1>, Dim<1>::Tensor>;
  template class ReplaceFieldList<Dim<1>, Dim<1>::SymTensor>;
                 
  template class ReplaceFieldList<Dim<2>, Dim<2>::Scalar>;
  template class ReplaceFieldList<Dim<2>, Dim<2>::Vector>;
  template class ReplaceFieldList<Dim<2>, Dim<2>::Vector3d>;
  template class ReplaceFieldList<Dim<2>, Dim<2>::Tensor>;
  template class ReplaceFieldList<Dim<2>, Dim<2>::SymTensor>;
                 
  template class ReplaceFieldList<Dim<3>, Dim<3>::Scalar>;
  template class ReplaceFieldList<Dim<3>, Dim<3>::Vector>;
  template class ReplaceFieldList<Dim<3>, Dim<3>::Vector3d>;
  template class ReplaceFieldList<Dim<3>, Dim<3>::Tensor>;
  template class ReplaceFieldList<Dim<3>, Dim<3>::SymTensor>;
}
