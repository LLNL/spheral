text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "PhysicsEvolvingMaterialLibrary.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class PhysicsEvolvingMaterialLibrary<Dim< %(ndim)s > >;
}
"""
