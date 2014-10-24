//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "IncrementFieldList.cc"

namespace Spheral {
  using FieldSpace::Field;
  template class IncrementFieldList<Dim<1>, Dim<1>::Scalar>;
  template class IncrementFieldList<Dim<1>, Dim<1>::Vector>;
  template class IncrementFieldList<Dim<1>, Dim<1>::Vector3d>;
  template class IncrementFieldList<Dim<1>, Dim<1>::Tensor>;
  template class IncrementFieldList<Dim<1>, Dim<1>::SymTensor>;
              
  template class IncrementFieldList<Dim<2>, Dim<2>::Scalar>;
  template class IncrementFieldList<Dim<2>, Dim<2>::Vector>;
  template class IncrementFieldList<Dim<2>, Dim<2>::Vector3d>;
  template class IncrementFieldList<Dim<2>, Dim<2>::Tensor>;
  template class IncrementFieldList<Dim<2>, Dim<2>::SymTensor>;
              
  template class IncrementFieldList<Dim<3>, Dim<3>::Scalar>;
  template class IncrementFieldList<Dim<3>, Dim<3>::Vector>;
  template class IncrementFieldList<Dim<3>, Dim<3>::Vector3d>;
  template class IncrementFieldList<Dim<3>, Dim<3>::Tensor>;
  template class IncrementFieldList<Dim<3>, Dim<3>::SymTensor>;
}
