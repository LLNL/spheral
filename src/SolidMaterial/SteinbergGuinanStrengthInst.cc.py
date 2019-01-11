text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SteinbergGuinanStrength.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class SteinbergGuinanStrength<Dim< %(ndim)s > >;
}
"""
