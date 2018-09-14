text = """
//------------------------------------------------------------------------------
// Instantiations.
//------------------------------------------------------------------------------
#include "Damage/TensorDamageModel.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class PhysicsSpace::TensorDamageModel<Dim< %(ndim)s > >;
}
"""
