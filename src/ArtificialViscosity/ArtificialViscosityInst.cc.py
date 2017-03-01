text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "ArtificialViscosity.cc"

namespace Spheral {
  namespace ArtificialViscositySpace {
    template class ArtificialViscosity< Dim< %(ndim)s > >;
  }
}
"""
