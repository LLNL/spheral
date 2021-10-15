text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "GSPH/Limiters/OspreLimiter.cc"

namespace Spheral {
  template class OspreLimiter<Dim< %(ndim)s > >;
}
"""
