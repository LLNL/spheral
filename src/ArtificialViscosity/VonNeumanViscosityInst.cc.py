text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "ArtificialConduction/VonNeumanViscosity.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace ArtificialViscositySpace {
    template class VonNeumanViscosity< Dim< %(ndim)s > >;
  }
}
"""
