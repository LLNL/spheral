text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "DataBase/CopyState.cc"

namespace Spheral {
  template class CopyState<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>;
  template class CopyState<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>;
  template class CopyState<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector3d>;
  template class CopyState<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>;
  template class CopyState<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>;
}
"""
