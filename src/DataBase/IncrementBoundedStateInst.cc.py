text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "DataBase/IncrementBoundedState.cc"

namespace Spheral {
  template class IncrementBoundedState<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>;
  template class IncrementBoundedState<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>;
  template class IncrementBoundedState<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector3d>;
  template class IncrementBoundedState<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>;
  template class IncrementBoundedState<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>;
  template class IncrementBoundedState<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor, Dim< %(ndim)s >::Scalar>;
}
"""
