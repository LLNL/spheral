text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "GSPH/Limiters/MinModLimiter.cc"

namespace Spheral {
  template class MinModLimiter<Dim< %(ndim)s > >;
}
"""
