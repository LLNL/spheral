text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "GSPH/Policies/PureReplaceWithStateFieldList.cc"

namespace Spheral {
  template class PureReplaceWithStateFieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>;
  template class PureReplaceWithStateFieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>;
  template class PureReplaceWithStateFieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector3d>;
  template class PureReplaceWithStateFieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>;
  template class PureReplaceWithStateFieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>;
}
"""
