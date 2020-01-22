text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "DataBase/ReplaceBoundedFieldList.cc"

namespace Spheral {
  template class ReplaceBoundedFieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>;
  template class ReplaceBoundedFieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>;
  template class ReplaceBoundedFieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector3d>;
  template class ReplaceBoundedFieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>;
  template class ReplaceBoundedFieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>;
  template class ReplaceBoundedFieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor, Dim< %(ndim)s >::Scalar>;
}
"""
