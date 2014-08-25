//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "IncrementState.cc"

namespace Spheral {
  using FieldSpace::Field;
  template class IncrementState<Dim<1>, Dim<1>::Scalar>;
  template class IncrementState<Dim<1>, Dim<1>::Vector>;
  template class IncrementState<Dim<1>, Dim<1>::Vector3d>;
  template class IncrementState<Dim<1>, Dim<1>::Tensor>;
  template class IncrementState<Dim<1>, Dim<1>::SymTensor>;
              
  template class IncrementState<Dim<2>, Dim<2>::Scalar>;
  template class IncrementState<Dim<2>, Dim<2>::Vector>;
  template class IncrementState<Dim<2>, Dim<2>::Vector3d>;
  template class IncrementState<Dim<2>, Dim<2>::Tensor>;
  template class IncrementState<Dim<2>, Dim<2>::SymTensor>;
              
  template class IncrementState<Dim<3>, Dim<3>::Scalar>;
  template class IncrementState<Dim<3>, Dim<3>::Vector>;
  template class IncrementState<Dim<3>, Dim<3>::Vector3d>;
  template class IncrementState<Dim<3>, Dim<3>::Tensor>;
  template class IncrementState<Dim<3>, Dim<3>::SymTensor>;
}
