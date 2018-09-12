text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "DataBase/IncrementFieldList.cc"

namespace Spheral {
  template class IncrementFieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>;
  template class IncrementFieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>;
  template class IncrementFieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector3d>;
  template class IncrementFieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>;
  template class IncrementFieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>;
}
"""
