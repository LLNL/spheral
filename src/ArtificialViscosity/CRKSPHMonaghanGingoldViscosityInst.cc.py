text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "ArtificialConduction/CRKSPHMonaghanGingoldViscosity.cc"

namespace Spheral {
  namespace ArtificialViscositySpace {
    template class CRKSPHMonaghanGingoldViscosity< Dim< %(ndim)s > >;
  }
}
"""
