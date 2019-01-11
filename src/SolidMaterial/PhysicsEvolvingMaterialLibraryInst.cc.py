text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SolidMaterial/PhysicsEvolvingMaterialLibrary.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class PhysicsEvolvingMaterialLibrary<Dim< %(ndim)s > >;
}
"""
