text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "GSPH/Limiters/VanLeerLimiter.cc"

namespace Spheral {
  template class VanLeerLimiter<Dim< %(ndim)s > >;
}
"""
