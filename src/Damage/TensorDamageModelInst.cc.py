text = """
//------------------------------------------------------------------------------
// Instantiations.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "TensorDamageModel.cc"

namespace Spheral {
  template class PhysicsSpace::TensorDamageModel<Dim< %(ndim)s > >;
}
"""
