text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "VonNeumanViscosity.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class VonNeumanViscosity< Dim< %(ndim)s > >;
}
"""
