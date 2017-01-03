text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "PhysicsEvolvingMaterialLibrary.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace SolidMaterial {
    template class PhysicsEvolvingMaterialLibrary<Dim< %(ndim)s > >;
  }
}
"""
