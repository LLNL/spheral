text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "DataBase/CopyFieldList.cc"

namespace Spheral {
  template class CopyFieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>;
  template class CopyFieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>;
  template class CopyFieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector3d>;
  template class CopyFieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>;
  template class CopyFieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>;
}
"""
