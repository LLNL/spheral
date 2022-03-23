text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "GSPH/Policies/PureReplaceFieldList.cc"

namespace Spheral {
  template class PureReplaceFieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>;
  template class PureReplaceFieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>;
  template class PureReplaceFieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector3d>;
  template class PureReplaceFieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>;
  template class PureReplaceFieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>;
}
"""
