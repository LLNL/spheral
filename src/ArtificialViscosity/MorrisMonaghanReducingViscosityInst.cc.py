text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "ArtificialViscosity/MorrisMonaghanReducingViscosity.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace ArtificialViscositySpace {
    template class MorrisMonaghanReducingViscosity< Dim< %(ndim)s > >;
  }
}
"""
