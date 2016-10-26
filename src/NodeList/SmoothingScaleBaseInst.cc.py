text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SmoothingScaleBase.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace NodeSpace {
    template class SmoothingScaleBase< Dim< %(ndim)s > >;
  }
}
"""
