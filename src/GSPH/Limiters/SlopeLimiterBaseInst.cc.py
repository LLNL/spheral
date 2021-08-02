text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "GSPH/Limiters/SlopeLimiterBase.cc"

namespace Spheral {
  template class SlopeLimiterBase<Dim< %(ndim)s > >;
}
"""
