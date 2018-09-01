text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "MorrisMonaghanReducingViscosity.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class MorrisMonaghanReducingViscosity< Dim< %(ndim)s > >;
}
"""
