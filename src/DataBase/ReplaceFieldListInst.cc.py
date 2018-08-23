text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "DataBase/ReplaceFieldList.cc"

namespace Spheral {
  template class ReplaceFieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>;
  template class ReplaceFieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>;
  template class ReplaceFieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector3d>;
  template class ReplaceFieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>;
  template class ReplaceFieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>;
}
"""
