text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "MorrisMonaghanReducingViscosity.cc"
#include "ArtificialConduction/Geometry/Dimension.hh"

namespace Spheral {
  namespace ArtificialViscositySpace {
    template class MorrisMonaghanReducingViscosity< Dim< %(ndim)s > >;
  }
}
"""
