text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "GSPH/Policies/IncrementSpecificFromTotalPolicy.cc"

namespace Spheral {
  template class IncrementSpecificFromTotalPolicy<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>;
  template class IncrementSpecificFromTotalPolicy<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>;
}
"""
