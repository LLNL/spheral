text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "GSPH/Limiters/VanAlbaLimiter.cc"

namespace Spheral {
  template class VanAlbaLimiter<Dim< %(ndim)s > >;
}
"""
