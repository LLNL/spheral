text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "GSPH/Limiters/SuperbeeLimiter.cc"

namespace Spheral {
  template class SuperbeeLimiter<Dim< %(ndim)s > >;
}
"""
