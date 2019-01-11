text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "CRKSPHMonaghanGingoldViscosity.cc"

namespace Spheral {
  template class CRKSPHMonaghanGingoldViscosity< Dim< %(ndim)s > >;
}
"""
