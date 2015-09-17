text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "CRKSPHMonaghanGingoldViscosity.cc"

namespace Spheral {
  namespace ArtificialViscositySpace {
    template class CRKSPHMonaghanGingoldViscosity< Dim< %(ndim)s > >;
  }
}
"""
