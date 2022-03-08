text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "GSPH/ReplaceWithRatioPolicy.cc"

namespace Spheral {
  template class ReplaceWithRatioPolicy<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>;
}
"""
