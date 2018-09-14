text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "DataBase/IncrementBoundedFieldList.cc"

namespace Spheral {
  template class IncrementBoundedFieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>;
  template class IncrementBoundedFieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>;
  template class IncrementBoundedFieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector3d>;
  template class IncrementBoundedFieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>;
  template class IncrementBoundedFieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>;
  template class IncrementBoundedFieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor, Dim< %(ndim)s >::Scalar>;
}
"""
