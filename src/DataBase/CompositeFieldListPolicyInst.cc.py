text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "DataBase/CompositeFieldListPolicy.cc"

namespace Spheral {
  template class CompositeFieldListPolicy<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>;
  template class CompositeFieldListPolicy<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>;
  template class CompositeFieldListPolicy<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector3d>;
  template class CompositeFieldListPolicy<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>;
  template class CompositeFieldListPolicy<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>;
}
"""
