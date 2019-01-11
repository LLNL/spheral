text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "FiniteVolumeViscosity.cc"

namespace Spheral {
  template class FiniteVolumeViscosity< Dim< %(ndim)s > >;
}
"""
