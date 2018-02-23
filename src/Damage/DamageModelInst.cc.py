text = """
//------------------------------------------------------------------------------
// Instantiations.
//------------------------------------------------------------------------------
#include "DamageModel.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class PhysicsSpace::DamageModel<Dim< %(ndim)s > >;
}
"""
