text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "DEM/ReplaceAndIncrementPairFieldList.cc"

namespace Spheral {
  template class ReplaceAndIncrementPairFieldList<Dim< %(ndim)s >,std::vector<Dim< %(ndim)s>::Scalar>>;
  template class ReplaceAndIncrementPairFieldList<Dim< %(ndim)s >,std::vector<Dim< %(ndim)s>::Vector>>;
}
"""
