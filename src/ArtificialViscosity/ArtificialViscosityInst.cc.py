text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "ArtificialViscosity/ArtificialViscosity.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace ArtificialViscositySpace {
    template class ArtificialViscosity< Dim< %(ndim)s > >;
  }
}
"""
