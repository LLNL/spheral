text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "ArtificialViscosity/VonNeumanViscosity.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class VonNeumanViscosity< Dim< %(ndim)s > >;
}
"""
