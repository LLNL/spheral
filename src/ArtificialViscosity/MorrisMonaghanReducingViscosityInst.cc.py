text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "ArtificialViscosity/MorrisMonaghanReducingViscosity.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class MorrisMonaghanReducingViscosity< Dim< %(ndim)s > >;
}
"""
