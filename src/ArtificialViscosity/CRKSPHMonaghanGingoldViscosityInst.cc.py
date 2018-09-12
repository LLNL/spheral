text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "ArtificialViscosity/CRKSPHMonaghanGingoldViscosity.cc"

namespace Spheral {
  template class CRKSPHMonaghanGingoldViscosity< Dim< %(ndim)s > >;
}
"""
