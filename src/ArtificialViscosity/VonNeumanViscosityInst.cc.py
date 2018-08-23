text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "ArtificialViscosity/VonNeumanViscosity.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace ArtificialViscositySpace {
    template class VonNeumanViscosity< Dim< %(ndim)s > >;
  }
}
"""
