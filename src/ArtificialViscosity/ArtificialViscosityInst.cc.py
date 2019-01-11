text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "ArtificialViscosity.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class ArtificialViscosity< Dim< %(ndim)s > >;
}
"""
