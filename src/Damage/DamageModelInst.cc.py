text = """
//------------------------------------------------------------------------------
// Instantiations.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "DamageModel.cc"

namespace Spheral {
  template class PhysicsSpace::DamageModel<Dim< %(ndim)s > >;
}
"""
