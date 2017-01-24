text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "VonNeumanViscosity.cc"

namespace Spheral {
  namespace ArtificialViscositySpace {
    template class VonNeumanViscosity< Dim< %(ndim)s > >;
  }
}
"""
