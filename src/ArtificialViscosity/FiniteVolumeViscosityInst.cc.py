text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "ArtificialViscosity/FiniteVolumeViscosity.cc"

namespace Spheral {
  template class FiniteVolumeViscosity< Dim< %(ndim)s > >;
}
"""
