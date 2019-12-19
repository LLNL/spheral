text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "DataBase/ReplaceBoundedState.cc"

namespace Spheral {
  template class ReplaceBoundedState<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>;
  template class ReplaceBoundedState<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>;
  template class ReplaceBoundedState<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector3d>;
  template class ReplaceBoundedState<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>;
  template class ReplaceBoundedState<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>;
  template class ReplaceBoundedState<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor, Dim< %(ndim)s >::Scalar>;
}
"""
