text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "MorrisMonaghanReducingViscosity.cc"

namespace Spheral {
  namespace ArtificialViscositySpace {
    template class MorrisMonaghanReducingViscosity< Dim< %(ndim)s > >;
  }
}
"""
