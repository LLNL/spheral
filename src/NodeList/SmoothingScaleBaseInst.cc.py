text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "NodeList/SmoothingScaleBase.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class SmoothingScaleBase< Dim< %(ndim)s > >;
}
"""
