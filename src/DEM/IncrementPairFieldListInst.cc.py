text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "DEM/IncrementPairFieldList.cc"

namespace Spheral {
  template class IncrementPairFieldList<Dim< %(ndim)s >, std::vector<int>>;
  template class IncrementPairFieldList<Dim< %(ndim)s >, std::vector<Dim< %(ndim)s >::Scalar>>;
  template class IncrementPairFieldList<Dim< %(ndim)s >, std::vector<Dim< %(ndim)s >::Vector>>;
}
"""
