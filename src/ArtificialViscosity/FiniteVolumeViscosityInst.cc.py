text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "FiniteVolumeViscosity.cc"

namespace Spheral {
  namespace ArtificialViscositySpace {
    template class FiniteVolumeViscosity< Dim< %(ndim)s > >;
  }
}
"""
