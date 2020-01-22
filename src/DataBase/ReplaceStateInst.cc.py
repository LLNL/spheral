text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "DataBase/ReplaceState.cc"

namespace Spheral {
  template class ReplaceState<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>;
  template class ReplaceState<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>;
  template class ReplaceState<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector3d>;
  template class ReplaceState<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>;
  template class ReplaceState<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>;
}
"""
