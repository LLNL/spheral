text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "GSPH/Policies/PureReplaceState.cc"

namespace Spheral {
  template class PureReplaceState<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>;
  template class PureReplaceState<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>;
  template class PureReplaceState<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector3d>;
  template class PureReplaceState<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>;
  template class PureReplaceState<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>;
}
"""
