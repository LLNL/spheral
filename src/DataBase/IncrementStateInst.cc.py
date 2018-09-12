text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "DataBase/IncrementState.cc"

namespace Spheral {
  template class IncrementState<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>;
  template class IncrementState<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>;
  template class IncrementState<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector3d>;
  template class IncrementState<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>;
  template class IncrementState<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>;
}
"""
