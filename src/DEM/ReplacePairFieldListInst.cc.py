text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "DEM/ReplacePairFieldList.cc"

namespace Spheral {
  template class ReplacePairFieldList<Dim< %(ndim)s >, std::vector<int>>;
  template class ReplacePairFieldList<Dim< %(ndim)s >, std::vector<Dim< %(ndim)s >::Scalar>>;
  template class ReplacePairFieldList<Dim< %(ndim)s >, std::vector<Dim< %(ndim)s >::Vector>>;
}
"""
