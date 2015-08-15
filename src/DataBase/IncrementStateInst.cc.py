text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "IncrementState.cc"

namespace Spheral {
  using FieldSpace::Field;
  template class IncrementState<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>;
  template class IncrementState<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>;
  template class IncrementState<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector3d>;
  template class IncrementState<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>;
  template class IncrementState<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>;
}
"""
