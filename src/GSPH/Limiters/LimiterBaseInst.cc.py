text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "GSPH/Limiters/LimiterBase.cc"

namespace Spheral {
  template class LimiterBase<Dim< %(ndim)s > >;
}
"""
