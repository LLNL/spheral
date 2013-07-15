//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "ReplaceBoundedState.cc"

namespace Spheral {
  template class ReplaceBoundedState<Dim<1>, Dim<1>::Scalar>;
  template class ReplaceBoundedState<Dim<1>, Dim<1>::Vector>;
  template class ReplaceBoundedState<Dim<1>, Dim<1>::Vector3d>;
  template class ReplaceBoundedState<Dim<1>, Dim<1>::Tensor>;
  template class ReplaceBoundedState<Dim<1>, Dim<1>::SymTensor>;
  template class ReplaceBoundedState<Dim<1>, Dim<1>::SymTensor, Dim<1>::Scalar>;
                 
  template class ReplaceBoundedState<Dim<2>, Dim<2>::Scalar>;
  template class ReplaceBoundedState<Dim<2>, Dim<2>::Vector>;
  template class ReplaceBoundedState<Dim<2>, Dim<2>::Vector3d>;
  template class ReplaceBoundedState<Dim<2>, Dim<2>::Tensor>;
  template class ReplaceBoundedState<Dim<2>, Dim<2>::SymTensor>;
  template class ReplaceBoundedState<Dim<2>, Dim<2>::SymTensor, Dim<2>::Scalar>;
                 
  template class ReplaceBoundedState<Dim<3>, Dim<3>::Scalar>;
  template class ReplaceBoundedState<Dim<3>, Dim<3>::Vector>;
  template class ReplaceBoundedState<Dim<3>, Dim<3>::Vector3d>;
  template class ReplaceBoundedState<Dim<3>, Dim<3>::Tensor>;
  template class ReplaceBoundedState<Dim<3>, Dim<3>::SymTensor>;
  template class ReplaceBoundedState<Dim<3>, Dim<3>::SymTensor, Dim<3>::Scalar>;
}
