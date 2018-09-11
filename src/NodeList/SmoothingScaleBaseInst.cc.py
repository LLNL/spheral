text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "NodeList/SmoothingScaleBase.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace NodeSpace {
    template class SmoothingScaleBase< Dim< %(ndim)s > >;
  }
}
"""
